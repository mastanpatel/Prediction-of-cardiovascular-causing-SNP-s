load("m3_lab2_save.RData")
library(plyr)

GWASout <- read.table("GWAA.txt", header = TRUE, colClasses = c("character", 
    rep("numeric", 4)))

GWASout$Neg_logP <- -log10(GWASout$p.value)

GWASout <- merge(GWASout, genoBim[, c("SNP", "chr", "position")])

GWASout <- arrange(GWASout, -Neg_logP)

GWASout$type <- "typed"

head(GWASout)
##          SNP   Estimate  Std.Error   t.value      p.value Neg_logP chr
## 1  rs1532625  0.2076611 0.03760084  5.522779 4.041152e-08 7.393495  16
## 2  rs3803768 -0.3128523 0.06769878 -4.621241 4.195818e-06 5.377183  17
## 3 rs16951746 -0.3219660 0.07856772 -4.097943 4.429065e-05 4.353688  15
## 4 rs17643302  0.2422581 0.06183232  3.917984 9.398109e-05 4.026960  17
## 5 rs17778044  0.1882179 0.04804567  3.917479 9.423620e-05 4.025782  17
## 6  rs4149504 -0.1824303 0.04664973 -3.910640 9.685440e-05 4.013881  16
##   position  type
## 1 57005301 typed
## 2 80872028 typed
## 3 68605486 typed
## 4 64770822 typed
## 5 68310164 typed
## 6 75565845 typed


GWASout$type <- "typed"

GWAScomb<-rbind.fill(GWASout, imputeOut)

head(GWAScomb)
##          SNP   Estimate  Std.Error   t.value      p.value Neg_logP chr
## 1  rs1532625  0.2076611 0.03760084  5.522779 4.041152e-08 7.393495  16
## 2  rs3803768 -0.3128523 0.06769878 -4.621241 4.195818e-06 5.377183  17
## 3 rs16951746 -0.3219660 0.07856772 -4.097943 4.429065e-05 4.353688  15
## 4 rs17643302  0.2422581 0.06183232  3.917984 9.398109e-05 4.026960  17
## 5 rs17778044  0.1882179 0.04804567  3.917479 9.423620e-05 4.025782  17
## 6  rs4149504 -0.1824303 0.04664973 -3.910640 9.685440e-05 4.013881  16
##   position  type
## 1 57005301 typed
## 2 80872028 typed
## 3 68605486 typed
## 4 64770822 typed
## 5 68310164 typed
## 6 75565845 typed
tail(GWAScomb)
##               SNP Estimate Std.Error t.value   p.value     Neg_logP chr
## 190467   rs857176       NA        NA      NA 0.9999661 1.472859e-05  16
## 190468   rs857177       NA        NA      NA 0.9999661 1.472859e-05  16
## 190469   rs860730       NA        NA      NA 0.9999661 1.472859e-05  16
## 190470   rs863228       NA        NA      NA 0.9999661 1.472859e-05  16
## 190471  rs9989429       NA        NA      NA 0.9999661 1.472859e-05  16
## 190472 rs72774978       NA        NA      NA 0.9999699 1.307962e-05  16
##        position    type
## 190467 13930841 imputed
## 190468 13930825 imputed
## 190469 13930975 imputed
## 190470 13931061 imputed
## 190471 13924679 imputed
## 190472  7661950 imputed

CETP <- map2gene("CETP", coords = genes, SNPs = GWAScomb)

CETP <- CETP[,c("SNP","p.value","Neg_logP","chr","position","type","gene")]

print(CETP)
##                SNP      p.value  Neg_logP chr position    type gene
## 1        rs1532625 4.041152e-08 7.3934948  16 57005301   typed CETP
## 16        rs289742 3.496236e-04 3.4563993  16 57017762   typed CETP
## 130       rs289715 4.047466e-03 2.3928168  16 57008508   typed CETP
## 32369    rs1800775 3.970058e-08 7.4012031  16 56995236 imputed CETP
## 32370    rs3816117 3.970058e-08 7.4012031  16 56996158 imputed CETP
## 32371    rs1532624 4.763363e-08 7.3220863  16 57005479 imputed CETP
## 32372    rs7205804 4.763363e-08 7.3220863  16 57004889 imputed CETP
## 32386   rs17231506 3.464271e-05 4.4603881  16 56994528 imputed CETP
## 32388     rs183130 3.464271e-05 4.4603881  16 56991363 imputed CETP
## 32391    rs3764261 3.464271e-05 4.4603881  16 56993324 imputed CETP
## 32392     rs821840 3.464271e-05 4.4603881  16 56993886 imputed CETP
## 32394     rs820299 4.474057e-05 4.3492985  16 57000284 imputed CETP
## 32420   rs12149545 9.991848e-05 4.0003542  16 56993161 imputed CETP
## 32445   rs12447839 2.061185e-04 3.6858831  16 56993935 imputed CETP
## 32446   rs12447924 2.061185e-04 3.6858831  16 56994192 imputed CETP
## 32448    rs4783962 2.061185e-04 3.6858831  16 56995038 imputed CETP
## 32480   rs34620476 3.516690e-04 3.4538660  16 56996649 imputed CETP
## 32481     rs708272 3.516690e-04 3.4538660  16 56996288 imputed CETP
## 32482     rs711752 3.516690e-04 3.4538660  16 56996211 imputed CETP
## 32525   rs12447620 3.573152e-04 3.4469486  16 57014319 imputed CETP
## 32526     rs158480 3.573152e-04 3.4469486  16 57008227 imputed CETP
## 32527     rs158617 3.573152e-04 3.4469486  16 57008287 imputed CETP
## 32754   rs11508026 1.932329e-03 2.7139189  16 56999328 imputed CETP
## 32755   rs12444012 1.932329e-03 2.7139189  16 57001438 imputed CETP
## 32756   rs12720926 1.932329e-03 2.7139189  16 56998918 imputed CETP
## 32757    rs4784741 1.932329e-03 2.7139189  16 57001216 imputed CETP
## 33224  rs112039804 4.055013e-03 2.3920078  16 57018856 imputed CETP
## 33225   rs12708985 4.055013e-03 2.3920078  16 57014610 imputed CETP
## 33226     rs736274 4.055013e-03 2.3920078  16 57009769 imputed CETP
## 37683     rs289746 3.208886e-02 1.4936456  16 57020205 imputed CETP
## 73171   rs12934552 2.569000e-01 0.5902359  16 57021433 imputed CETP
## 84766   rs12708983 3.234698e-01 0.4901663  16 57014411 imputed CETP
## 116938  rs12598522 5.233758e-01 0.2811864  16 57022352 imputed CETP
## 116939  rs56315364 5.233758e-01 0.2811864  16 57021524 imputed CETP

write.csv(CETP, "CETP.csv", row.names = FALSE)  # save for future use

save(imputed, target, rules, phenodata, support, genes, impCETP, impCETPgeno, 
    imputeOut, GWAA, map2gene, GWASout, GWAScomb, CETP, file = "m4_lab1_save.RData")