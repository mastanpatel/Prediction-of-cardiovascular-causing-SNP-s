load("lab3_save.RData")
source("http://bioconductor.org/biocLite.R")
biocLite("snpStats")
## 
## The downloaded binary packages are in
##  /var/folders/_9/8f93yrbn70v__gsn6550wbqr0000gn/T//Rtmp2scphX/downloaded_packages
biocLite("gdsfmt")
## 
## The downloaded binary packages are in
##  /var/folders/_9/8f93yrbn70v__gsn6550wbqr0000gn/T//Rtmp2scphX/downloaded_packages
library(snpStats)
library(gdsfmt)
hardy <- 10^-6     
CADcontrols <- as.character(clinical[ clinical$CAD==0, 'FamID' ])
snpsum.colCont <- col.summary(genotype[CADcontrols,])
HWEuse <- with(snpsum.colCont, !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(hardy/2) ) ) )
rm(snpsum.colCont)
HWEuse[is.na(HWEuse)] <- FALSE          
cat(ncol(genotype)-sum(HWEuse),"SNPs will be removed due to high HWE.\n") 
## 765 SNPs will be removed due to high HWE.
genotype <- genotype[,HWEuse]
print(genotype)                           
## A SnpMatrix with  1400 rows and  381257 columns
## Row names:  10002 ... 11596 
## Col names:  rs4579145 ... rs946221

save(genotype, genoBim, clinical, file= "lab4_save.RData")