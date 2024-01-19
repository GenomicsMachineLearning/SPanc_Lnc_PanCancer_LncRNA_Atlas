setwd("~/Downloads/Acral")
# a1, a2, a5, a6, c1 --> v2
a1<-Read10X_h5("GSM6622292_AM1_filtered_feature_bc_matrix.h5")
a2<-Read10X_h5("GSM6622293_AM2_filtered_feature_bc_matrix.h5")
a3<-Read10X_h5("GSM6622294_AM3_pre_filtered_feature_bc_matrix.h5")
a3_post<-Read10X_h5("GSM6622295_AM3_post_filtered_feature_bc_matrix.h5")
a4<-Read10X_h5("GSM6622296_AM4_filtered_feature_bc_matrix.h5")
c1<-Read10X_h5("GSM6622299_CM1_filtered_feature_bc_matrix.h5")
c1_lym<-Read10X_h5("GSM6622302_CM1_lym_filtered_feature_bc_matrix.h5")
c2<-Read10X_h5("GSM6622300_CM2_filtered_feature_bc_matrix.h5")

a1_s<-CreateSeuratObject(a1,project = "Acral1")
a2_s<-CreateSeuratObject(a2,project = "Acral2")
a3_s<-CreateSeuratObject(a3,project = "Acral3_Pre")
a3_post_s<-CreateSeuratObject(a3_post,project = "Acral3_Post")
a4_s<-CreateSeuratObject(a4,project = "Acral4")
c1_s<-CreateSeuratObject(c1,project = "Cutaneous1")
c1_l_s<-CreateSeuratObject(c1_lym,project = "Cutaneous1_Lymph")
c2_s<-CreateSeuratObject(c2,project = "Cutaneous2")


sample.list=c("a1_s","a2_s","a3_s","a3_post_s","a4_s","c1_s","c1_l_s","c2_s")
sample.list <- list(A1=a1_s,A2=a2_s,A3_Pre=a3_s,A3_Post=a3_post_s,A4=a4_s,Cutaneous1=c1_s,Cutaneous1_Lymph=c1_l_s,Cutaneous2=c2_s)
sample.list <- lapply(X = sample.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = sample.list)
anchors <- FindIntegrationAnchors(object.list = sample.list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
pdf("integrated.pdf", width = 12, height = 8.5)
p1 + p2
dev.off()
saveRDS(combined,"seurat_obj.rds")

mel<-c("MITF,PMEL,TYR","DCT","MLANA","PMEL","APOC1","S100A1") #0,15,19  16,11 1,2,5,6
mel_cycle<-c("UBE2C","NUSAP1","MKI67","CENPF") #7,18
endo<-c("VWF","PECAM1") #10
fibro<-c("COL3A1","COL1A2","COL1A1","LUM") #14,13
CD8T<-c("CD8A","HAVCR2","LAG3","PD1","TIGIT","CTLA4","HOPX") #3. 23 ,4
mono<-c("CD14","LYZ","CD74","CD68","CD79A") #17,22,21
B<-c("MS4A1") #12
CD4T<-c("CD4","FOXP3","IL2")
NK<-c("GNLY","FGFBP2","FCGR3A","KLRD1","KLRF1")
NKT<-c("CD3E","CD3D","GZMB","XCL2","IFNG","CCL4","NKG7","GZMA","GZMK")
macrophages<-c('CD68','CD163','CD14','CD11b','CD206','CD80','CD86','CD16','CD64','CCL18','CD115','CD11c','CD32','HLA-DR','MRC1','MSR1','GCA','Pf4')

Mel_markers<-c("PIK3CD","PIK3R3","NRAS","AKT3","RAF1","MITF","PIK3CB","PIK3CA","FGF12","PDGFRA","FGF5","EGF","FGF2","PDGFC","FGF10","PIK3R1","FGF1","PDGFRB","FGF18","E2F3","CDKN1A","PDGFA","EGFR","HGF","CDK6","PIK3CG","MET","BRAF","FGF20","FGF17","FGFR1","CDKN2A","PTEN","FGF8","HRAS","BAD","CCND1","FGF19","FGF4","FGF3","PDGFD","FGF23","FGF6","CDK4","MDM2","IGF1","FGF9","RB1","FGF14","AKT1","FGF7","MAP2K1","IGF1R","MAPK3","CDH1","FGF11","TP53","PIK3R5","FGF22","MAP2K2","PIK3R2","AKT2","FGF21","E2F1","MAPK1","PDGFB","ARAF","FGF16","FGF13")
cancer<-c("ALK","AFP","BCL2","B2M","BTA","BRCA1","BRCA2","BRAF","CD117","CEA","CD19","CD22","FLT3","MYC","KRAS","JAK2","IRF4","IDH1","FGFR2","P53","RAS")
classical_monocytes<-c('CD14','CD11b','CD68','HLA-DR','CD33','CD11c','CD123','CD15','CD3D','CD3E','CD3G','CD3Z','CD66b','VCAN','S100A12','CXCL8','S100A8','S100A9','LYZ','CST3','Elane1')
nonclassical_monocytes<-c('CD14','CD16','CD11b','CD68','HLA-DR','CD33','CD11c','CD123','CD15','CD3D','CD3E','CD3G','CD3Z','CD66b','FCGR3A','CDKN1C','LST1','FCER1G','MS4A7','RHOC','S100A8','S100A9','CST3','C1QC','Elane1')
intermediate_monocytes<-c('CD14','CD16','CD11b','CD68','HLA-DR','CD33','CD11c','CD123','CD15','CD3D','CD3E','CD3G','CD3Z','CD66b','IL1B','S100A8','S100A9','CST3','C1QC','Elane1')

DotPlot(combined,features=c(mel,mel_cycle,endo,fibro,CD8T,B,mono,NK,CD4T)) + RotatedAxis()
DotPlot(combined,features=c(Mel_markers,cancer)) + RotatedAxis()


combined <- RenameIdents(object = combined, `0` = "Melanoma", `1` = "Melanoma", `2` = "Melanoma", `5` = "Melanoma", `6` = "Melanoma", `11` = "Melanoma", `15` = "Melanoma",`16` = "Melanoma", `19` = "Melanoma",`18` = "Cycling Melanoma",`7` = "Cycling Melanoma",`10` = "Endothelial cells",`13` = "Fibroblasts",`14` = "Fibroblasts", `3` = "Memory CD8+ T cells",`23` = "CD8 T cells",`4` = "CD8+ NK T-like cells",`21` = "Monocytes",`22` = "Monocytes",`17` = "Monocytes",`12` = "B cells")
DimPlot(combined, label = TRUE)

## automatic annotation scType
source("sctype_score_.R")
source("gene_sets_prepare.R")
source("auto_detect_tissue_type.R")

all.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoMultiBarHeatmap(combined, assay = 'RNA', features = top_10$gene, group.by='Idents', additional.group.by = 'orig.ident')

#In addition, provide a tissue type your data belongs to:
# DB file
db_ = "ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)
es.max = sctype_score(scRNAseqData = combined[["integrated"]]@scale.data, scaled = TRUE,gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(combined@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(combined@meta.data[combined@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(combined@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"


## OVERLAYING IDENTIFIED CELL-TYPES ON PLOTS
combined@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  combined@meta.data$customclassif[combined@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")

# UMAP and spatial plot with inferred cell-types
#Idents(combined)=combined@meta.data$customclassif
DimPlot(combined, reduction = "umap", label = TRUE, cols = ccolss, group.by = "customclassif")


#auto-detect tissue
#####################################################################################################
# read combined_seurat.rds
subset<-readRDS("combined_seurat.rds")
subset<-subset(subset,idents = c("Acral1","Acral2","Acral3_Post","Acral3_Pre","Acral4","Cut1","Cut1_Lym","Cut2" ))
subset<-NormalizeData(subset)
subset<-ScaleData(subset)


subset <- FindVariableFeatures(subset, selection.method = "vst", nfeatures = 2000)
subset <- RunPCA(object = subset, reduction.name="pca",features = VariableFeatures(subset))
subset <- FindNeighbors(object=subset, reduction="pca", graph.name="gene_snn",dims=1:30, )
subset <- FindClusters(object=subset,graph.name="gene_snn",resolution=0.2)
subset <- RunUMAP(subset,reduction.name="umap",reduction="pca",dims=1:30,check_duplicates = FALSE)
DimPlot(subset, reduction = "umap")

pdf("umap.pdf")
DimPlot(subset, reduction = "umap")
dev.off()
pdf("sample_umap.pdf")
DimPlot(subset, group.by = "ident")
dev.off()
