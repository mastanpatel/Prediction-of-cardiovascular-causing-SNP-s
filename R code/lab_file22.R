load("m2_lab1_save.RData")
library(snpStats)
library(downloader)
info <- download("https://www.mtholyoke.edu/courses/afoulkes/Data/statsTeachR/chr16_1000g_CEU.info", 
    destfile = "chr16_1000g_CEU.info")
ped <- download("https://www.mtholyoke.edu/courses/afoulkes/Data/statsTeachR/chr16_1000g_CEU.ped", 
    destfile = "chr16_1000g_CEU.ped")
thougeno <- read.pedfile("chr16_1000g_CEU.ped", snps = "chr16_1000g_CEU.info", 
    which = 1)
genoMatrix <- thougeno$genotypes
support <- thougeno$map
colnames(support) <- c("SNP", "position", "A1", "A2")
head(support)
##           SNP position A1 A2
## 1 rs140769322    60180  3  2
## 2 rs188810967    60288  2  1
## 3  rs76368850    60291  2  4
## 4 rs185537431    60778  3  1
## 5 rs542544747    60842  2  1
## 6   rs4021615    61349  1  3

presSnps <- colnames(genotype)
presSnps <- colnames(genotype)
presDatChr <- genoBim[genoBim$SNP %in% presSnps & genoBim$chr == 16, ]
targetSnps <- presDatChr$SNP
is.present <- colnames(genoMatrix) %in% targetSnps
missing <- genoMatrix[, !is.present]
print(missing)
## A SnpMatrix with  99 rows and  386404 columns
## Row names:  CEU_1 ... CEU_99 
## Col names:  rs140769322 ... rs111706106
present <- genoMatrix[, is.present]
print(present)
## A SnpMatrix with  99 rows and  12047 columns
## Row names:  CEU_1 ... CEU_99 
## Col names:  rs41340949 ... rs7196459

pos.pres <- support$position[is.present]
pos.miss <- support$position[!is.present]
rules <- snp.imputation(present, missing, pos.pres, pos.miss)
## SNPs tagged by a single SNP: 68562
## SNPs tagged by multiple tag haplotypes (saturated model): 137495
rules <- rules[can.impute(rules)]
cat("Imputation rules for", length(rules), "SNPs were estimated\n") 
## Imputation rules for 206057 SNPs were estimated

# Quality control for imputation certainty and MAF
r2threshold <- 0.7
minor <- 0.01
rules <- rules[imputation.r2(rules) >= r2threshold]
cat(length(rules),"imputation rules remain after uncertain imputations were removed\n")  
## 159085 imputation rules remain after uncertain imputations were removed
rules <- rules[imputation.maf(rules) >= minor]
cat(length(rules),"imputation rules remain after MAF filtering\n") 
## 159085 imputation rules remain after MAF filtering
# Obtain posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules, target, as.numeric=FALSE)
print(imputed)  
## A SnpMatrix with  1400 rows and  159085 columns
## Row names:  10002 ... 11596 
## Col names:  rs28436661 ... rs62053708

rm(genoMatrix)
rm(missing)
rm(present)
save(genotype, genoBim, clinical, pcs, imputed, target, rules, support, file="m2_lab2_save.RData")