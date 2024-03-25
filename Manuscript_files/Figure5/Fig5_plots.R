####################################### Manhattan plot #################################
m2<-read.csv("~/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Vis15_HN//C1/Manhattan/qqman_data.txt",sep="\t",header=TRUE)
uTARs<- grep("^chr", m2$gene)
uTARs <- m2[uTARs, ]
cuTARs <- subset(uTARs, pval > 0.05 | pval < -0.05)
uTARs<-uTARs$gene
manhattan1(m2, chr = "Chrom", bp = "start", p = "Log2FC", snp = "gene",col = c("#e5e4e2", "#bebebe"),highlight1 = uTARs, highlight2 = cuTARs$gene, logp = FALSE,ylim=c(-8,8), ylab="average log2FC",suggestiveline = F, genomewideline = F)
abline(h = 2.5, col = "black", lty = 3)  # Red line
abline(h = -2.5, col = "black", lty = 3)
