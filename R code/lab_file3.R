load("lab2_save.RData")
library(snpStats)
source("http://bioconductor.org/biocLite.R")
biocLite("SNPRelate")
## 
## The downloaded binary packages are in
##  /var/folders/_9/8f93yrbn70v__gsn6550wbqr0000gn/T//Rtmp8kW974/downloaded_packages
library(SNPRelate)                     
library(plyr)
snpsum.row <- row.summary(genotype)
MAF <- snpsum.col$MAF
callmatrix <- !is.na(genotype)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(snpsum.row, Heterozygosity*(ncol(genotype))*Call.rate)
snpsum.row$hetF <- 1-(hetObs/hetExp)
head(snpsum.row)
##       Call.rate Certain.calls Heterozygosity         hetF
## 10002 0.9826528             1      0.3281344 -0.022019308
## 10004 0.9891708             1      0.3231777 -0.006829093
## 10005 0.9918277             1      0.3239852 -0.008817640
## 10007 0.9860113             1      0.3241203 -0.009882026
## 10008 0.9824172             1      0.3231505 -0.008550116
## 10009 0.9912570             1      0.3205117  0.002331055
hetcutoff <- 0.1   
sampleuse <- with(snpsum.row, abs(hetF) <= hetcutoff)
sampleuse[is.na(sampleuse)] <- FALSE  
cat(nrow(genotype)-sum(sampleuse)
## 0 subjects will be removed due to low sample call rate or inbreeding coefficient.

genotype <- genotype[sampleuse,]
clinical<- clinical[rownames(genotype), ]

print(MAF[13])
## [1] 0.02607143

sampcall <- 0.95   
sampleuse1 <- with(snpsum.row, !is.na(Call.rate) & Call.rate > sampcall)
sampleuse1[is.na(sampleuse1)] <- FALSE   
cat(nrow(genotype)-sum(sampleuse1), "subjects will be removed due to low sample call rate or inbreeding coefficient.\n")
## 0 subjects will be removed due to low sample call rate or inbreeding coefficient.
genotype <- genotype[sampleuse1,]
clinical<- clinical[ rownames(genotype), ]ld.thresh <- 0.2    
kin.thresh <- 0.1   
snpgdsBED2GDS("GWAS_data.bed", "GWAS_data.fam", "GWAS_data.bim","GWAS_data.gds")
## Start snpgdsBED2GDS ...
##  BED file: "GWAS_data.bed" in the SNP-major mode (Sample X SNP)
##  FAM file: "GWAS_data.fam", DONE.
##  BIM file: "GWAS_data.bim", DONE.
## Thu Aug 27 11:58:23 2015     store sample id, snp id, position, and chromosome.
##  start writing: 1401 samples, 500000 SNPs ...
##      Thu Aug 27 11:58:23 2015    0%
##      Thu Aug 27 11:58:32 2015    100%
## Thu Aug 27 11:58:33 2015     Done.
## Optimize the access efficiency ...
## Clean up the fragments of GDS file:
##  open the file "GWAS_data.gds" (size: 179677809).
##  # of fragments in total: 39.
##  save it to "GWAS_data.gds.tmp".
##  rename "GWAS_data.gds.tmp" (size: 179677557).
##  # of fragments in total: 18.
genofile <- openfn.gds("GWAS_data.gds", readonly = FALSE)
gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
add.gdsn(genofile, "sample.id", gds.ids, replace = TRUE)

set.seed(1000)
geno.sample.ids <- rownames(genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
                          sample.id = geno.sample.ids, 
                          snp.id = colnames(genotype)) 
## SNP pruning based on LD:
## Excluding 117978 SNPs on non-autosomes
## Excluding 0 SNP (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
## Working space: 1401 samples, 382022 SNPs
##  Using 1 (CPU) core
##  Sliding window: 500000 basepairs, Inf SNPs
##  |LD| threshold: 0.2
## Chromosome 4: 75.65%, 24341/32174
## Chromosome 6: 77.35%, 24251/31353
## Chromosome 9: 77.01%, 18363/23846
## Chromosome 18: 75.73%, 11477/15156
## Chromosome 3: 76.32%, 26881/35222
## Chromosome 17: 75.83%, 8783/11583
## Chromosome 19: 76.03%, 5071/6670
## Chromosome 11: 76.12%, 19474/25583
## Chromosome 10: 76.53%, 21411/27978
## Chromosome 1: 75.08%, 31002/41294
## Chromosome 13: 76.03%, 15115/19880
## Chromosome 8: 76.89%, 21622/28122
## Chromosome 12: 76.43%, 18631/24377
## Chromosome 22: 74.48%, 4875/6545
## Chromosome 2: 75.94%, 32479/42770
## Chromosome 16: 75.62%, 12164/16085
## Chromosome 5: 77.26%, 25166/32572
## Chromosome 7: 77.56%, 20985/27058
## Chromosome 15: 76.05%, 11429/15028
## Chromosome 14: 76.46%, 12447/16280
## Chromosome 20: 76.10%, 10076/13241
## Chromosome 21: 76.75%, 5513/7183
## 381556 SNPs are selected in total.
snpset.ibd <- unlist(snpSUB, use.names=FALSE)
cat(length(snpset.ibd),"will be used in IBD analysis\n")
## 381556 will be used in IBD analysis

ibd <- snpgdsIBDMoM(genofile, kinship=TRUE,
                    sample.id = geno.sample.ids,
                    snp.id = snpset.ibd,
                    num.thread = 1)
## IBD analysis (PLINK method of moment) on SNP genotypes:
## Excluding 118444 SNPs on non-autosomes
## Excluding 0 SNP (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
## Working space: 1401 samples, 381556 SNPs
##  Using 1 (CPU) core
## PLINK IBD:   the sum of all working genotypes (0, 1 and 2) = 251556755
## PLINK IBD:   Thu Aug 27 11:58:45 2015    0%
## PLINK IBD:   Thu Aug 27 11:59:19 2015    18%
## PLINK IBD:   Thu Aug 27 11:59:52 2015    36%
## PLINK IBD:   Thu Aug 27 12:00:25 2015    54%
## PLINK IBD:   Thu Aug 27 12:00:59 2015    72%
## PLINK IBD:   Thu Aug 27 12:01:31 2015    91%
## PLINK IBD:   Thu Aug 27 12:01:47 2015    100%
ibdcoeff <- snpgdsIBDSelection(ibd)    
head(ibdcoeff)
##     ID1   ID2        k0         k1     kinship
## 1 10002 10004 0.9651755 0.03482454 0.008706134
## 2 10002 10005 0.9694550 0.03054498 0.007636244
## 3 10002 10007 0.9327868 0.06721323 0.016803307
## 4 10002 10008 0.9218918 0.07810819 0.019527049
## 5 10002 10009 0.9666561 0.03334394 0.008335986
## 6 10002 10010 0.9715478 0.02845221 0.007113054

ibdcoeff <- ibdcoeff[ ibdcoeff$kinship >= kin.thresh, ]
related.samples <- NULL
while ( nrow(ibdcoeff) > 0 ) {
    sample.counts <- arrange(count(c(ibdcoeff$ID1, ibdcoeff$ID2)), -freq)
    rm.sample <- sample.counts[1, 'x']
    cat("Removing sample", as.character(rm.sample), 'too closely related to', 
        sample.counts[1, 'freq'],'other samples.\n')

    ibdcoeff <- ibdcoeff[ibdcoeff$ID1 != rm.sample & ibdcoeff$ID2 != rm.sample,]
    related.samples <- c(as.character(rm.sample), related.samples)
}
## Removing sample 10347 too closely related to 1 other samples.

genotype <- genotype[ !(rownames(genotype) %in% related.samples), ]
clinical <- clinical[ !(clinical$FamID %in% related.samples), ]

geno.sample.ids <- rownames(genotype)

cat(length(related.samples), "similar samples removed due to correlation coefficient >=", 
    kin.thresh,"\n") 
## 1 similar samples removed due to correlation coefficient >= 0.1
print(genotype)        
## A SnpMatrix with  1400 rows and  382022 columns
## Row names:  10002 ... 11596 
## Col names:  rs4579145 ... rs946221

pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids,  snp.id = snpset.ibd, num.thread=1)
## Principal Component Analysis (PCA) on SNP genotypes:
## Excluding 118444 SNPs on non-autosomes
## Excluding 0 SNP (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
## Working space: 1400 samples, 381556 SNPs
##  Using 1 (CPU) core
## PCA: the sum of all working genotypes (0, 1 and 2) = 251378768
## PCA: Thu Aug 27 12:02:05 2015    0%
## PCA: Thu Aug 27 12:02:35 2015    16%
## PCA: Thu Aug 27 12:03:05 2015    32%
## PCA: Thu Aug 27 12:03:35 2015    48%
## PCA: Thu Aug 27 12:04:05 2015    64%
## PCA: Thu Aug 27 12:04:35 2015    80%
## PCA: Thu Aug 27 12:05:05 2015    96%
## PCA: Thu Aug 27 12:05:11 2015    100%
## PCA: Thu Aug 27 12:05:11 2015    Begin (eigenvalues and eigenvectors)
## PCA: Thu Aug 27 12:05:12 2015    End (eigenvalues and eigenvectors)
pctab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[,1],    
                    PC2 = pca$eigenvect[,2],    
                    stringsAsFactors = FALSE)
plot(pctab$PC2, pctab$PC1, xlab="Principal Component 2", ylab="Principal Component 1", 
     main = "Ancestry Plot")

	 closefn.gds(genofile)
save(genotype, snpsum.col, snpsum.row, genofile, genoBim, clinical, file= "lab3_save.RData")