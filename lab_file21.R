load("lab4_save.RData")
library(SNPRelate)                    
genofile <- snpgdsOpen("GWAS_data.gds", readonly = TRUE)ld.thresh <- 0.2
set.seed(1000)
geno.sample.ids <- rownames(genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh, sample.id = geno.sample.ids, 
    snp.id = colnames(genotype))
## SNP pruning based on LD:
## Excluding 118743 SNPs on non-autosomes
## Excluding 0 SNP (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
## Working space: 1400 samples, 381257 SNPs
##  Using 1 (CPU) core
##  Sliding window: 500000 basepairs, Inf SNPs
##  |LD| threshold: 0.2
## Chromosome 4: 75.51%, 24294/32174
## Chromosome 6: 77.23%, 24214/31353
## Chromosome 9: 76.84%, 18323/23846
## Chromosome 18: 75.59%, 11456/15156
## Chromosome 3: 76.15%, 26823/35222
## Chromosome 17: 75.66%, 8764/11583
## Chromosome 19: 75.86%, 5060/6670
## Chromosome 11: 75.98%, 19438/25583
## Chromosome 10: 76.33%, 21356/27978
## Chromosome 1: 74.92%, 30937/41294
## Chromosome 13: 75.85%, 15079/19880
## Chromosome 8: 76.75%, 21583/28122
## Chromosome 12: 76.26%, 18591/24377
## Chromosome 22: 74.36%, 4867/6545
## Chromosome 2: 75.80%, 32421/42770
## Chromosome 16: 75.51%, 12146/16085
## Chromosome 5: 77.04%, 25093/32572
## Chromosome 7: 77.42%, 20949/27058
## Chromosome 15: 75.90%, 11407/15028
## Chromosome 14: 76.31%, 12424/16280
## Chromosome 20: 75.99%, 10062/13241
## Chromosome 21: 76.64%, 5505/7183
## 380792 SNPs are selected in total.
snpset.pca <- unlist(snpSUB, use.names = FALSE)
cat(length(snpset.pca), "\n")
## 380792

pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids,  snp.id = snpset.pca, num.thread=1)
## Principal Component Analysis (PCA) on SNP genotypes:
## Excluding 119208 SNPs on non-autosomes
## Excluding 0 SNP (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
## Working space: 1400 samples, 380792 SNPs
##  Using 1 (CPU) core
## PCA: the sum of all working genotypes (0, 1 and 2) = 251238746
## PCA: Thu Aug 27 12:59:58 2015    0%
## PCA: Thu Aug 27 13:00:28 2015    14%
## PCA: Thu Aug 27 13:00:58 2015    29%
## PCA: Thu Aug 27 13:01:29 2015    44%
## PCA: Thu Aug 27 13:01:59 2015    58%
## PCA: Thu Aug 27 13:02:29 2015    73%
## PCA: Thu Aug 27 13:02:59 2015    88%
## PCA: Thu Aug 27 13:03:23 2015    100%
## PCA: Thu Aug 27 13:03:23 2015    Begin (eigenvalues and eigenvectors)
## PCA: Thu Aug 27 13:03:24 2015    End (eigenvalues and eigenvectors)
pcs <- data.frame(FamID = pca$sample.id, pca$eigenvect[,1 : 10],
                  stringsAsFactors = FALSE)
colnames(pcs)[2:11]<-paste("pc", 1:10, sep = "")
print(head(pcs))
##   FamID         pc1          pc2          pc3          pc4          pc5
## 1 10002 -0.01121786  0.015282391  0.012335017  0.017869885 -0.013561539
## 2 10004  0.01263937  0.007136085 -0.009394939 -0.029605670  0.017640584
## 3 10005  0.01656767 -0.007822581  0.034279450 -0.008364165  0.043822928
## 4 10007  0.01179249 -0.001340544  0.014388948  0.002448683  0.008521085
## 5 10008  0.01587414  0.003683702  0.013407115 -0.002069099  0.009863999
## 6 10009  0.01342526  0.002620955  0.005132100 -0.044515366  0.043797908
##            pc6          pc7          pc8          pc9         pc10
## 1  0.015099899  0.005363108 -0.007476769 -0.003362287 -0.007697605
## 2  0.062301372 -0.005341143  0.004102591  0.040899908  0.003655463
## 3 -0.002068763  0.015205957 -0.011775232  0.005064646  0.003447803
## 4  0.017414424  0.015268401  0.037495537  0.016946863  0.010266512
## 5  0.049090388  0.007638991  0.010654394 -0.013350077  0.007254919
## 6  0.012536302  0.045885623 -0.033450390 -0.045956099  0.009781019

closefn.gds(genofile)
save(genotype, genoBim, clinical, pcs, file="m2_lab1_save.RData")