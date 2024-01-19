# Main figures
########################### TCGA HM #################################
#setwd("~/Downloads/PhD/Analysis/ST/QUAN/RNAscope/new/")
T<-read.csv("~/Downloads/Manuscript_final_figs/TCGA_HM_renamed_data.txt", sep="\t", header=TRUE)
MT <- T[, -1]
rownames(MT)<-T$Sample.ID
Tanno<-as.data.frame(read.csv("~/Downloads/Manuscript_final_figs/TCGA_HM_Metadata_renamed.txt", sep="\t", header=TRUE))
Tannom<- as.data.frame(as.matrix(Tanno[, -1]))
rownames(Tannom)<-Tanno$uTAR
pheatmap(log(MT+0.0001),cluster_rows=FALSE, annotation_row = Tannom, angle=45)


### Supplementary 
#### nanopore umi-bc downsamplin
d<-read.csv("UMII_BC_sampling",sep="\t",header=TRUE)
df <- d %>% filter(No..of.UMIs < 1000)
ggplot(df, aes(x = as.factor(No..of.UMIs), y = No.of.uTARs)) + geom_boxplot()+  theme_minimal() +labs(x = "No. of UMIs per spot", y = "No.of uTARs")

