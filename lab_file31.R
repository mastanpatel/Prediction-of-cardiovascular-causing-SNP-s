load("m2_lab2_save.RData")  
library(snpStats)
library(plyr)
library(GenABEL)
library(doParallel)

genodata <- genotype
#Print the number of SNPs to be checked
cat(paste(ncol(genodata), " SNPs included in analysis.\n"))
## 381257  SNPs included in analysis.
#create text file for GWAA output to be written to
columns<-c("SNP", "Estimate", "Std.Error", "t-value", "p-value")
write.table(t(columns), "GWAA.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

phenoSub <- merge(clinical,pcs)     
phenoSub$phenotype <- rntransform(as.numeric(phenoSub$hdl), family="gaussian")
phenoSub <- rename(phenoSub, replace=c(FamID="id"))
phenoSub<-phenoSub[!is.na(phenoSub$phenotype),]
phenodata <- phenoSub
genodata <- genodata[as.character(phenodata$id),]
cat(nrow(genodata), "samples included in analysis.\n")
## 1308 samples included in analysis.
phenodata$hdl <- NULL
phenodata$ldl <- NULL
phenodata$tg <- NULL
phenodata$CAD <- NULL
print(head(phenodata))
##      id sex age        pc1          pc2          pc3          pc4
## 2 10004   2  50 0.01263937  0.007136085 -0.009394939 -0.029605670
## 3 10005   1  55 0.01656767 -0.007822581  0.034279450 -0.008364165
## 4 10007   1  52 0.01179249 -0.001340544  0.014388948  0.002448683
## 5 10008   1  58 0.01587414  0.003683702  0.013407115 -0.002069099
## 6 10009   1  59 0.01342526  0.002620955  0.005132100 -0.044515366
## 7 10010   1  54 0.01627605  0.011860496  0.018480130 -0.014176319
##            pc5          pc6          pc7          pc8          pc9
## 2  0.017640584  0.062301372 -0.005341143  0.004102591  0.040899908
## 3  0.043822928 -0.002068763  0.015205957 -0.011775232  0.005064646
## 4  0.008521085  0.017414424  0.015268401  0.037495537  0.016946863
## 5  0.009863999  0.049090388  0.007638991  0.010654394 -0.013350077
## 6  0.043797908  0.012536302  0.045885623 -0.033450390 -0.045956099
## 7 -0.003986151 -0.094898367  0.004855911 -0.014761558  0.018721010
##           pc10  phenotype
## 2  0.003655463 -2.2874211
## 3  0.003447803 -0.4742508
## 4  0.010266512  0.8850182
## 5  0.007254919 -0.1636196
## 6  0.009781019  0.3952471
## 7 -0.019222291  1.7105947

par(mfrow=c(1,2))
hist(as.numeric(phenoSub$hdl), main=NULL, xlab="HDL")
hist(phenodata$phenotype, main=NULL, xlab="Transformed HDL")

GWAA <- function(genodata, phenodata, filename = NULL, append = FALSE, workers = getOption("mc.cores", 
    2L), flip = TRUE, select.snps = NULL, hosts = NULL) {
    
    if (is.null(hosts)) {
        cl <- makeCluster(workers)
    } else {
        cl <- makeCluster(hosts, "PSOCK")
    }
    show(cl)
    registerDoParallel(cl)
    
    # Function that will change which allele is counted (major or minor)
    flip.matrix <- function(x) {
        zero2 <- which(x == 0)
        two0 <- which(x == 2)
        x[zero2] <- 2
        x[two0] <- 0
        return(x)
    }
    
    foreach(part = 1:nSplits) %do% {
        genoNum <- as(genodata[, snp.start[part]:snp.stop[part]], "numeric")
        # flip.matrix function employed
        if (isTRUE(flip)) 
            genoNum <- flip.matrix(genoNum)
        rsVec <- colnames(genoNum)
        res <- foreach(snp.name = rsVec, .combine = "rbind") %dopar% {
            a <- summary(glm(phenotype ~ . - id, family = gaussian, data = cbind(phenodata, 
                snp = genoNum[, snp.name])))
            a$coefficients["snp", ]
        }
        
        write.table(cbind(rsVec, res), filename, append = TRUE, quote = FALSE, 
            col.names = FALSE, row.names = FALSE)
        cat(sprintf("GWAS SNPs %s-%s (%s%% finished)\n", snp.start[part], snp.stop[part], 
            100 * part/nSplits))
    }
    
    stopCluster(cl)
    return(print("Done."))
}


SNPs_sub <- genoBim$SNP[genoBim$chr==15 | genoBim$chr==16 
                        | genoBim$chr==17]
genodata_sub <- genodata[,colnames(genodata)%in%SNPs_sub]

nSNPs <- ncol(genodata_sub)
nSplits <- 10
genosplit <- ceiling(nSNPs/nSplits) # number of SNPs in each subset
snp.start <- seq(1, nSNPs, genosplit) # index of first SNP in group
snp.stop <- pmin(snp.start+genosplit-1, nSNPs) # index of last SNP in group

start <- Sys.time()
GWAA(genodata_sub, phenodata, filename="GWAA.txt")
## socket cluster with 2 nodes on host 'localhost'
## GWAS SNPs 1-3237 (10% finished)
## GWAS SNPs 3238-6474 (20% finished)
## GWAS SNPs 6475-9711 (30% finished)
## GWAS SNPs 9712-12948 (40% finished)
## GWAS SNPs 12949-16185 (50% finished)
## GWAS SNPs 16186-19422 (60% finished)
## GWAS SNPs 19423-22659 (70% finished)
## GWAS SNPs 22660-25896 (80% finished)
## GWAS SNPs 25897-29133 (90% finished)
## GWAS SNPs 29134-32368 (100% finished)
## [1] "Done."
end <- Sys.time()
print(end-start)
## Time difference of 2.71885 mins

save(genoBim, imputed, target, rules, phenodata, genodata_sub, support, GWAA, 
    file = "m3_lab1_save.RData")