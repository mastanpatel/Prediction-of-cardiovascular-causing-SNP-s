load("m3_lab1_save.RData")
library(snpStats)
library(plyr)
rownames(phenodata) <- as.character(phenodata$id)

imp <- snp.rhs.tests(phenotype ~ sex + age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + 
    pc7 + pc8 + pc9 + pc10, family = "Gaussian", data = phenodata, snp.data = target, 
    rules = rules)

results <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), stringsAsFactors = FALSE)
results <- results[!is.na(results$p.value), ]

write.csv(results, "impute.csv", row.names = FALSE)

imputeOut <- merge(results, support[, c("SNP", "position")])
imputeOut$chr <- 16

imputeOut$type <- "imputed"

imputeOut$Neg_logP <- -log10(imputeOut$p.value)

imputeOut <- arrange(imputeOut, p.value)
print(head(imputeOut))
##          SNP      p.value position chr    type Neg_logP
## 1  rs1800775 3.970058e-08 56995236  16 imputed 7.401203
## 2  rs3816117 3.970058e-08 56996158  16 imputed 7.401203
## 3  rs1532624 4.763363e-08 57005479  16 imputed 7.322086
## 4  rs7205804 4.763363e-08 57004889  16 imputed 7.322086
## 5 rs12933833 2.112374e-05 56697684  16 imputed 4.675229
## 6 rs11076159 2.400306e-05 56670827  16 imputed 4.619733

map2gene <- function(gene, coords, SNPs, extend.boundary = 5000) {
  coordsSub <- coords[coords$gene == gene,] 

  coordsSub$start <- coordsSub$start - extend.boundary 
  coordsSub$stop <- coordsSub$stop + extend.boundary

  SNPsub <- SNPs[SNPs$position >= coordsSub$start & SNPs$position <= coordsSub$stop &
                 SNPs$chr == coordsSub$chr,] 

  return(data.frame(SNPsub, gene = gene, stringsAsFactors = FALSE))
}

# Read in file containing protein coding genes coords
library(downloader)
ProdCod <- download("https://www.mtholyoke.edu/courses/afoulkes/Data/statsTeachR/ProCodgene_coords.csv", 
    destfile = "ProCodgene_coords.csv")
genes <- read.csv("ProCodgene_coords.csv", stringsAsFactors = FALSE)

# Subset for CETP SNPs
impCETP <- map2gene("CETP", coords = genes, SNPs = imputeOut)

# Filter for CETP SNP genotypes
impCETPgeno <- imputed[, impCETP$SNP]

save(genoBim, imputed, target, rules, phenodata, support, genes, impCETP, impCETPgeno, 
    imputeOut, GWAA, map2gene, file = "m3_lab2_save.RData")