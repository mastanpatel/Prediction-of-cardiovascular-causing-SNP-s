source("http://bioconductor.org/biocLite.R")
biocLite("snpStats")
## 
## The downloaded binary packages are in
##  /var/folders/q5/brz1dtz54w99g1mft67shggm0000gp/T//Rtmp2RngR5/downloaded_packages
library(snpStats)
library(RCurl)


# read in data in two steps:

library(downloader)
bedFile<- download("https://www.mtholyoke.edu/courses/afoulkes/Data/statsTeachR/GWAS_data.bed",
destfile = "GWAS_data.bed")
bimFile<- download("https://www.mtholyoke.edu/courses/afoulkes/Data/statsTeachR/GWAS_data.bim", destfile = "GWAS_data.bim")
famFile<- download("https://www.mtholyoke.edu/courses/afoulkes/Data/statsTeachR/GWAS_data.fam",
destfile = "GWAS_data.fam")

# step 2: read in the local files using read.plink()
data <- read.plink("GWAS_data.bed","GWAS_data.bim","GWAS_data.fam",na.strings=("-9"))

# what class is "data"
class(data)
## [1] "list"
# how many elements are in "data"
length(data)
## [1] 3
bim <- data$map
head(bim)
##            chromosome   snp.name cM  position allele.1 allele.2
## rs4579145           4  rs4579145  0  78303781        T        C
## rs2768995           6  rs2768995  0   6911315        C        T
## rs10125738          9 rs10125738  0 101066136        A        G
## rs888263           18   rs888263  0  12750253        G        A
## rs7639361           3  rs7639361  0  80138205        A        C
## rs2430512          17  rs2430512  0  69675940        T        C
(bed <- data$genotypes)
## A SnpMatrix with  1401 rows and  500000 columns
## Row names:  10002 ... 11596 
## Col names:  rs4579145 ... rs946221
fam <- data$fam
head(fam)
##       pedigree member father mother sex affected
## 10002    10002      1      0      0   1        1
## 10004    10004      1      0      0   2        1
## 10005    10005      1      0      0   1        2
## 10007    10007      1      0      0   1        1
## 10008    10008      1      0      0   1        2
## 10009    10009      1      0      0   1        2
library(RCurl)
clinicalURL <- getURL("https://www.mtholyoke.edu/courses/afoulkes/Data/statsTeachR/GWAS_clinical.csv")
clinical <- read.csv(text = clinicalURL, colClasses = c("character", rep("factor",
2), rep("numeric", 2)))
rownames(clinical) <- clinical$FamID
print(head(clinical))
##       FamID CAD sex age  tg  hdl  ldl
## 10002 10002   1   1  60  NA <NA> <NA>
## 10004 10004   1   2  50  55   23   75
## 10005 10005   1   1  55 105   37   69
## 10007 10007   1   1  52 314   54  108
## 10008 10008   1   1  58 161   40   94
## 10009 10009   1   1  59 171   46   92

genotype <- data$genotypes
genoBim <- data$map
colnames(genoBim)<-c("chr", "SNP", "gen.dist", "position", "A1", "A2")

save(genotype,genoBim,clinical,file="lab1_save.RData")