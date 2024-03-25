####################################### Manhattan plot #################################
m2<-read.csv("/Users/ritanya/Downloads/Acral/v2/wdgeR_pseudo/qqman_data.txt",sep="\t",header=TRUE)
uTARs<- grep("^uTAR", m2$gene)
uTARs <- m2[uTARs, ]
#cuTARs <- subset(uTARs, pval > 0.05 | pval < -0.05)
cuTARs<-c("uTAR97","uTAR117","uTAR650","uTAR915","uTAR1224","uTAR1249","uTAR1274","uTAR157","uTAR295","uTAR337")
#"cuTAR23054","cuTAR24389","cuTAR141997","cuTAR201736","cuTAR280980","cuTAR291293","cuTAR293520","cuTAR30875","cuTAR59116","cuTAR67350")
uTARs<-uTARs$gene
manhattan1(m2, chr = "chr", bp = "start", p = "avglogFC", snp = "gene",col = c("#e5e4e2", "#bebebe"),highlight1 = uTARs, highlight2 = cuTARs, logp = FALSE,ylim=c(-8,8), ylab="average log2FC",suggestiveline = F, genomewideline = F)
abline(h = 2.5, col = "black", lty = 3)  # Red line
abline(h = -2.5, col = "black", lty = 3)
