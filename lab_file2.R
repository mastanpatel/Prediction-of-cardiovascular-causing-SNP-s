load("lab1_save.RData")
library(snpStats)
## Loading required package: survival
## Loading required package: Matrix
## 
## Attaching package: 'Matrix'
## 
## The following objects are masked from 'package:base':
## 
##     crossprod, tcrossprod
snpsum.col <- col.summary(genotype)
head(snpsum.col)
##            Calls Call.rate Certain.calls       RAF         MAF
## rs4579145   1400 0.9992862             1 0.8235714 0.176428571
## rs2768995   1074 0.7665953             1 0.9790503 0.020949721
## rs10125738  1401 1.0000000             1 0.9768023 0.023197716
## rs888263    1382 0.9864383             1 0.6649783 0.335021708
## rs7639361   1396 0.9964311             1 0.6840974 0.315902579
## rs2430512   1385 0.9885796             1 0.9989170 0.001083032
##                    P.AA        P.AB      P.BB       z.HWE
## rs4579145  0.0285714286 0.295714286 0.6757143  0.65809530
## rs2768995  0.0037243948 0.034450652 0.9618250 -5.24953586
## rs10125738 0.0007137759 0.044967880 0.9543183 -0.29013170
## rs888263   0.1041968162 0.461649783 0.4341534  1.34207569
## rs7639361  0.0988538682 0.434097421 0.4670487  0.16261598
## rs2430512  0.0000000000 0.002166065 0.9978339  0.04034939

call <- 0.95
use <- with(snpsum.col, (!is.na(Call.rate) & Call.rate >= call))
use[is.na(use)] <- FALSE              
cat(ncol(genotype)-sum(use),"SNPs will be removed due to low call rate.\n") 
## 54583 SNPs will be removed due to low call rate.

genotype <- genotype[,use]
snpsum.col <- snpsum.col[use,]
print(genotype)         
## A SnpMatrix with  1401 rows and  445417 columns
## Row names:  10002 ... 11596 
## Col names:  rs4579145 ... rs946221

minor <- 0.01
use1 <- with(snpsum.col, (!is.na(MAF) & MAF > minor) )
use1[is.na(use1)] <- FALSE               
cat(ncol(genotype)-sum(use1),"SNPs will be removed due to low MAF .\n"  ) 
## 63395 SNPs will be removed due to low MAF .

genotype <- genotype[,use1]
snpsum.col <- snpsum.col[use1,]
print(genotype)
## A SnpMatrix with  1401 rows and  382022 columns
## Row names:  10002 ... 11596 
## Col names:  rs4579145 ... rs946221

save(genotype, snpsum.col,genoBim, clinical, file= "lab2_save.RData")