library(doParallel)
library(snpStats)
library(GenABEL)
load("m4_lab1_save.RData")

gwas1 <- GWASout[GWASout$SNP%in%colnames(target),]

genosub <- target[as.character(phenodata$id),]

phenosub <- phenodata[,c("id","phenotype")]
nSNPs <- ncol(genosub)
nSplits <- 10
genosplit <- ceiling(nSNPs/nSplits) # number of SNPs in each subset
snp.start <- seq(1, nSNPs, genosplit) # index of first SNP in group
snp.stop <- pmin(snp.start+genosplit-1, nSNPs) # index of last SNP in group

columns<-c("SNP", "Estimate", "Std.Error", "t-value", "p-value")
write.table(t(columns), "GWAA_unadj.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

gwas2 <- GWAA(genodata=genosub, phenodata=phenosub, filename="GWAA_unadj.txt")
## socket cluster with 2 nodes on host 'localhost'
## GWAS SNPs 1-1217 (10% finished)
## GWAS SNPs 1218-2434 (20% finished)
## GWAS SNPs 2435-3651 (30% finished)
## GWAS SNPs 3652-4868 (40% finished)
## GWAS SNPs 4869-6085 (50% finished)
## GWAS SNPs 6086-7302 (60% finished)
## GWAS SNPs 7303-8519 (70% finished)
## GWAS SNPs 8520-9736 (80% finished)
## GWAS SNPs 9737-10953 (90% finished)
## GWAS SNPs 10954-12164 (100% finished)
## [1] "Done."

# read in model results
GWASoutUnadj <- read.table("GWAA_unadj.txt", header = TRUE, colClasses = c("character", 
    rep("numeric", 4)))

# QQ plot for adjusted model
lambdaAdj <- estlambda(GWASout$t.value^2, plot = TRUE, method = "median", main = "Adjusted")

# QQ plot for unadjuested model
lambdaUnadj <- estlambda(GWASoutUnadj$t.value^2, plot = TRUE, method = "median", 
    main = "Unadjusted")
	
	c("Unadjusted model lambda:", lambdaUnadj$estimate)
## [1] "Unadjusted model lambda:" "1.05405027919361"
c("Adjusted model lambda:", lambdaAdj$estimate)
## [1] "Adjusted model lambda:" "1.03449915443322"

save(imputed, target, rules, phenodata, support, genes, impCETP, impCETPgeno, 
    imputeOut, GWAA, map2gene, GWASout, GWAScomb, CETP, file = "m4_lab2_save.RData")