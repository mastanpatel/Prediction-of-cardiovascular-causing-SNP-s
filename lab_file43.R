library(downloader)
library(LDheatmap)
library(rtracklayer)
library(ggplot2)
library(snpStats)
library(plyr)
library(devtools)
install_url("http://cran.r-project.org/src/contrib/Archive/postgwas/postgwas_1.11.tar.gz")
library(postgwas)
load("m4_lab2_save.RData")

code <- download("https://www.mtholyoke.edu/courses/afoulkes/Data/statsTeachR/GWAS_manhattan.R",
                           destfile="GWAS_manhattan.R")

source("GWAS_manhattan.R")

GWAS_Manhattan(GWAScomb)

subgen <- cbind(target[,colnames(target) %in% CETP$SNP], impCETPgeno)

# Subset SNPs for certain genotypes
certain <- apply(as(subgen, 'numeric'), 2, function(x) { all(x %in% c(0,1,2,NA)) })
subgen <- subgen[,certain]

CETP <- CETP[CETP$SNP %in% colnames(subgen),]

CETP <- arrange(CETP, position)

subgen <- subgen[ ,order(match(colnames(subgen),CETP$SNP))]

# Create LDheatmap
ld <- ld(subgen, subgen, stats="R.squared")

ll <- LDheatmap(ld, CETP$position, flip=TRUE, name="myLDgrob", title=NULL)

llplusgenes <- LDheatmap.addGenes(ll, chr = "chr16", genome = "hg19", genesLocation = 0.01)
plot.new()
llQplot2<-LDheatmap.addGrob(llplusgenes, rectGrob(gp = gpar(col = "white")),height = .34)
pushViewport(viewport(x = 0.483, y= 0.76, width = .91 ,height = .4))

grid.draw(ggplotGrob({
  qplot(position, Neg_logP, data = CETP, xlab="", ylab = "Negative Log P-value", 
        xlim = range(CETP$position), asp = 1/10, color = factor(type), 
        colour=c("#000000", "#D55E00")) + 
    theme(axis.text.x = element_blank(),
          axis.title.y = element_text(size = rel(0.75)), legend.position = "none", 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    scale_color_manual(values = c("red", "black"))
}))

snps<-data.frame(SNP=c("rs1532625"))

GWAScomb <- rename(GWAScomb, c(p.value="P", chr="CHR", position="BP"))

# Edit biomartConfigs so regionalplot function
# pulls from human genome build 37/hg19
myconfig <- biomartConfigs$hsapiens
myconfig$hsapiens$gene$host <- "grch37.ensembl.org"
myconfig$hsapiens$gene$mart <- "ENSEMBL_MART_ENSEMBL"
myconfig$hsapiens$snp$host <- "grch37.ensembl.org"
myconfig$hsapiens$snp$mart <- "ENSEMBL_MART_SNP"
Now we can use the following code to create a regional association plot, which will display the same data as the LD plot, but on a larger scale and in a slightly different format.

regionalplot(snps, GWAScomb, biomart.config = myconfig, 
             window.size = 400000, 
             draw.snpname = data.frame(
                snps = c("rs1532625"), 
                text = c("rs1532625"),
                angle = c(20),
                length = c(1), 
                cex = c(0.8)
                ),
ld.options = list(
  gts.source = 2, 
  max.snps.per.window = 2000, 
  rsquare.min = 0.8, 
  show.rsquare.text = FALSE
),
out.format = list(file = NULL, panels.per.page = 1))