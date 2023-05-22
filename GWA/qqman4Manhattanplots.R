
library('qqman')
library('tidyverse')
library(calibrate)  
  
  
Gap.dat <- read.csv("GAPIT.MLM.Color.qqBOI.csv", head=T)
For.Gap.dat <- Gap.dat  
 -log10(For.Gap.dat$P[1])
  
  
  Chr03 <- which(For.Gap.dat$CHR == 3)
  Chr03 <- For.Gap.dat[Chr03,]
  Ordf <- Chr03[which(Chr03$BP > 5000000 & Chr03$BP < 5200000),]
  Chr05 <- which(For.Gap.dat$CHR == 5)
  Chr05 <- For.Gap.dat[Chr05,]
  Ydf <- Chr05[which(Chr05$BP > 29935615 & Chr05$BP < 30065947),]   
  Chr07 <- which(For.Gap.dat$CHR == 7)
  Chr07 <- For.Gap.dat[Chr07,]
  Y2df <- Chr07[which(Chr07$BP > 38656562 & Chr07$BP < 39560030),]   
  OrSNPs <- as.character(Ordf$SNP)
  YSNPs <- as.character(Ydf$SNP)
  Y2SNPs <- as.character(Y2df$SNP)
  
  SNPsofFun <- c(OrSNPs, YSNPs, Y2SNPs)
  
  fix(manhattan)
  # set textxy 
  
  pdf(file="RG.Taproot.pdf", bg="white", width=10, height=6)
  manhattan(For.Gap.dat, main="Taproot Color", ylim=c(0,20),cex.lab=1.25, cex.axis = 1.25, cex.main=1.5,
            col=c("#C5050C","grey"), suggestiveline=F, genomewideline=7.86, chrlabs=c("1","2","3","4","5","6","7","8","9"),
            annotatePval  = 0.000000007, highlight = SNPsofFun)
  dev.off()
  
pdf(file="OB.Taproot.pdf", bg="white", width=10, height=6)
  manhattan(For.Gap.dat, main="Taproot Color", ylim=c(0,20),cex.lab=1.25, cex.axis = 1.25, cex.main=1.5,
            col=c("Orange","Blue"), suggestiveline=F, genomewideline=7.86, chrlabs=c("1","2","3","4","5","6","7","8","9"),
            annotatePval  = 0.000000007, highlight = SNPsofFun)
  dev.off()


Gap.dat <- read.csv("GAPIT.MLM.Carotene.qqBOI.csv", head=T)
For.Gap.dat <- Gap.dat  
-log10(For.Gap.dat$P[1])
  
Chr03 <- which(For.Gap.dat$CHR == 3)
Chr03 <- For.Gap.dat[Chr03,]
Ordf <- Chr03[which(Chr03$BP > 5000000 & Chr03$BP < 5200000),]
Chr05 <- which(For.Gap.dat$CHR == 5)
Chr05 <- For.Gap.dat[Chr05,]
Ydf <- Chr05[which(Chr05$BP > 29935615 & Chr05$BP < 30065947),]   
Chr07 <- which(For.Gap.dat$CHR == 7)
Chr07 <- For.Gap.dat[Chr07,]
Y2df <- Chr07[which(Chr07$BP > 38656562 & Chr07$BP < 39560030),]   
OrSNPs <- as.character(Ordf$SNP)
YSNPs <- as.character(Ydf$SNP)
Y2SNPs <- as.character(Y2df$SNP)

SNPsofFun <- c(OrSNPs, YSNPs, Y2SNPs)

  pdf(file="RG.alpabeta.pdf", bg="white", width=10, height=6)
  manhattan(For.Gap.dat, main="Percent ?? + ?? carotene", ylim=c(0,15),cex.lab=1.25, cex.axis = 1.25, cex.main=1.5,
            col=c("#C5050C","grey"), suggestiveline=F, genomewideline=7.86, chrlabs=c("1","2","3","4","5","6","7","8","9"),
            annotatePval  = 0.000000007, highlight = SNPsofFun)
  dev.off()
  
  pdf(file="OB.alpabeta.pdf", bg="white", width=10, height=6)
  manhattan(For.Gap.dat, main="Percent ?? + ?? carotene", ylim=c(0,15),cex.lab=1.25, cex.axis = 1.25, cex.main=1.5,
            col=c("Orange","Blue"), suggestiveline=F, genomewideline=7.86, chrlabs=c("1","2","3","4","5","6","7","8","9"),
            annotatePval  = 0.000000007, highlight = SNPsofFun)
  dev.off()  

  
  
  Gap.dat <- read.csv("GAPIT.MLM.Lutein.qqBOI.csv", head=T)
  For.Gap.dat <- Gap.dat  
  -log10(For.Gap.dat$P[1])
  
  Chr03 <- which(For.Gap.dat$CHR == 3)
  Chr03 <- For.Gap.dat[Chr03,]
  Ordf <- Chr03[which(Chr03$BP > 5000000 & Chr03$BP < 5200000),]
  Chr05 <- which(For.Gap.dat$CHR == 5)
  Chr05 <- For.Gap.dat[Chr05,]
  Ydf <- Chr05[which(Chr05$BP > 29935615 & Chr05$BP < 30065947),]   
  Chr07 <- which(For.Gap.dat$CHR == 7)
  Chr07 <- For.Gap.dat[Chr07,]
  Y2df <- Chr07[which(Chr07$BP > 38656562 & Chr07$BP < 39560030),]   
  OrSNPs <- as.character(Ordf$SNP)
  YSNPs <- as.character(Ydf$SNP)
  Y2SNPs <- as.character(Y2df$SNP)
  
  SNPsofFun <- c(OrSNPs, YSNPs, Y2SNPs)
  
  pdf(file="RG.Lutein.Gap.pdf", bg="white", width=10, height=6)
  manhattan(For.Gap.dat,main=" Percent Lutein", ylim=c(0,15),cex.lab=1.25, cex.axis = 1.25, cex.main=1.5,col=c("#C5050C","grey"), suggestiveline=F, genomewideline=7.86, chrlabs=c("1","2","3","4","5","6","7","8","9"),
            annotatePval  = 0.000000007, highlight = SNPsofFun)
  dev.off()
  
  pdf(file="OB.Lutein.pdf", bg="white", width=10, height=6)
  manhattan(For.Gap.dat, main="Percent Lutein", ylim=c(0,15),cex.lab=1.25, cex.axis = 1.25, cex.main=1.5,
            col=c("Orange","Blue"), suggestiveline=F, genomewideline=7.86, chrlabs=c("1","2","3","4","5","6","7","8","9"),
            annotatePval  = 0.000000007, highlight = SNPsofFun)
  dev.off()   