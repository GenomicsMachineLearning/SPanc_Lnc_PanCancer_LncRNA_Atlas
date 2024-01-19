library(edgeR)
so<-readRDS("v2_combined_seurat_obj.rds")
y<-Seurat2PB(so,sample="orig.ident",cluster="gene_cell_types")
summary(y$samples$lib.size)

keep.samples<-y$samples$lib.size > 50000 
table(keep.samples) 
#mel_samples<-c("Post Treatment_clusterMelanoma","Post Treatment_clusterCycling Melanoma","Pre Treatment_clusterCycling Melanoma","Pre Treatment_clusterMelanoma") 

y<-y[,keep.samples]
#Keep only uTARs for DE analysis. Skip to calculate FC for gene-co-exp
uTARs <- rownames(y)[grep("^uTAR", rownames(y))]
y<-y[uTARs,,keep=FALSE]

#normalization
y<-normLibSizes(y)

#MDS plots
cluster<-as.factor(y$samples$cluster) 
pdf("MDS_uTAR.pdf")
plotMDS(y,pch=16,col=c(2:8)[cluster],main="MDS") 
legend("topleft",legend=paste0(levels(cluster)), pch=14,col=2:8,cex=0.8)
dev.off()

#design matrix
donor<-factor(y$sample$sample)
design<-model.matrix(~cluster+donor)
colnames(design)<-gsub("donor","",colnames(design))
colnames(design)<-gsub("cluster","",colnames(design))
colnames(design)[1]<-"Int"

#dispersion estimates
#Dispersion Estimation
y<-estimateDisp(y,design,robust=TRUE)
y$common.dispersion #0.14  #0.2852. #V3 0.05g
pdf("dispersion.pdf")
plotBCV(y)
dev.off()

fit<-glmQLFit(y,design,robust=TRUE) 
pdf("QLdisp.pdf")
plotQLDisp(fit)
dev.off()


## Marker gene identification
ncls<-nlevels(cluster) 
contr<-rbind(matrix(1/(1-ncls),ncls,ncls), +matrix(0,ncol(design)-ncls,ncls)) 
diag(contr)<-1 
contr[1,]<-0 
rownames(contr)<-colnames(design) 
colnames(contr)<-levels(cluster)
contr

# Qasi-likelihood f-test for each comprison
qlf<-list() 
for(i in 1:ncls){ 
  qlf[[i]]<-glmQLFTest(fit,contrast=contr[,i]) 
  qlf[[i]]$comparison<-paste0(levels(donor)[i],"_vs_others") 
}

non_mel<-c("B cells","CD8+ NKT-like cells","Endothelial cells","Fibroblasts","Memory CD4+ T-cells","Monocytes","NK cells")        
mel<-c("Cycling Melanoma","Melanoma")

qlf[[11]]<-glmQLFTest(fit,contrast=contr[,mel]) 
qlf[[11]]$comparison<-paste0("Mel_vs_others") 
saveRDS(qlf,"qlf_uTARs.rds")
write.table(qlf[[11]]$table,"mel_vs_nonMel_all_genes_uTARs_qlf_table.txt",sep="\t",quote=FALSE,row.names=TRUE)

#No.of DE genes for each comparison
dt<-lapply(lapply(qlf,decideTestsDGE),summary) 
dt.all<-do.call("cbind",dt) 
head(dt.all)
write.table(dt.all,"cluster_wise_uTAR_stats.txt",sep="\t",quote=FALSE,row.names=TRUE)


# Top 20 DE uTARs per cluster
top<-10
topMarkers<-list()
for(i in 1:ncls){ 
  ord<-order(qlf[[i]]$table$PValue,decreasing=FALSE) 
  up<-qlf[[i]]$table$logFC>0 
  topMarkers[[i]]<-rownames(y)[ord[up][1:top]] 
} 
topMarkers<-unique(unlist(topMarkers)) 
write.table(topMarkers,"topmarkers_edger.txt",sep="\t",row.names = TRUE, quote = FALSE)


# Heat Map
lcpm<-cpm(y,log=TRUE) 
annot<-data.frame(cluster=cluster) 
rownames(annot)<-colnames(y) 
ann_colors<-list(cluster=2:11) 
names(ann_colors$cluster)<-levels(cluster)

colnames(lcpm)<-gsub("cluster","",colnames(lcpm))
pdf("uTAR_HM.pdf")
pheatmap::pheatmap(lcpm[topMarkers,],breaks=seq(-2,2,length.out=101),color=colorRampPalette(c("blue","white","red"))(100),scale="row",
                   cluster_cols=TRUE,border_color="NA",fontsize_row=5,treeheight_row=70,treeheight_col=70,cutree_cols=7,
                   clustering_method="ward.D2",show_colnames=FALSE,
                   annotation_col=annot,annotation_colors=ann_colors)
dev.off()


o<-read.csv("order_to_keep_pheatmap.txt",sep="\t",header=FALSE)
desired_order<-o$V1
lcpm_o <- lcpm[,match(desired_order, colnames(lcpm)), drop = FALSE]


rownames(annot)<-gsub("cluster","",rownames(annot))
annot_o <- annot[match(desired_order, rownames(annot)), , drop = FALSE]




pdf("uTAR_HM_ordered_label.pdf")
pheatmap::pheatmap(lcpm_o[topMarkers,],breaks=seq(-2,2,length.out=101),color=colorRampPalette(c("blue","white","red"))(100),scale="row",
                   cluster_cols=FALSE,border_color="NA",fontsize_row=5,
                   show_colnames=FALSE,
                   annotation_col=annot_o,annotation_colors=ann_colors)
#grid.text(levels(annot_o$cluster), x = c(0.25, 0.6), y = c(0.89, 0.89),gp = gpar(fontsize = 10, rot = 45))
dev.off()

Error in hclust(d, method = method) : 
  NA/NaN/Inf in foreign function call (arg 10)
Calls: <Anonymous> -> cluster_mat -> hclust
Execution halted





pdf("feature_plots_top_celltype_uTARs.pdf",height=10,width=8)
FeaturePlot(so, features = c("uTAR650","uTAR550","uTAR817","uTAR504","uTAR127","uTAR952"), max.cutoff = 3,
            cols = c("#f2f2f2", "red"),order=TRUE, spl)
dev.off()

genes_to_plot<-c("uTAR650","uTAR550","uTAR817","uTAR504","uTAR127","uTAR952")
plots <- lapply(genes_to_plot, function(gene) {
  FeaturePlot(so, features = gene, pt.size = 0.5)
})
combined_plot <- grid.arrange(grobs = plots, ncol = 2)  # Adjust ncol as needed
combined_plot <- combined_plot + theme(legend.position = "bottom") +
  guides(fill = FALSE)


pdf("feature_plots_top_celltype_uTARs_so.pdf",height=8,width=8)
FeaturePlot(so, features = c("uTAR1155","uTAR728","uTAR588","uTAR1043"), max.cutoff = 3,
            cols = c("#f2f2f2", "red"),order=TRUE)
dev.off()


Aall.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1, features = uTARs)
Atop10 <- Aall.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("DE_uTARs_seurat_HM.pdf",width=12,height=8)
DoHeatmap(combined, features = Atop10$gene)
dev.off()
