library(data.table)
library(Seurat)
bed<-read.csv("~/Downloads/all_cuTARs_refined_after_gencode47_clean.bed",sep="\t",header=FALSE)
#noDirMat38_HC<-("~/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Vis15_HN/C1/HNC_uTARs_final_mat.txt")
noDirMat38_HC<-("~/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Vis15_HN/C1/HNC_uTARs_final_mat.txt")
Sample_HMM_HC=fread(noDirMat38_HC,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
rownames(Sample_HMM_HC)<-tolower(rownames(Sample_HMM_HC))
Sample_HMM_HC <- as.data.frame(Sample_HMM_HC) # force to data frame
rownames(Sample_HMM_HC) <- Sample_HMM_HC$GENE # make the name of rows GENEs
Sample_HMM_HC<-Sample_HMM_HC[rownames(Sample_HMM_HC) %in% bed$V4,]
Sample_HMM_HC <- Sample_HMM_HC[,-1] # take out first column
uTAR_seurat_HC<-CreateSeuratObject(counts = Sample_HMM_HC, project = "Vis15_HeadNeck_C",min.cells = 3)


HC<- Read10X(data.dir ="~/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Vis15_HN/C1/filtered_feature_bc_matrix/", gene.column = 2, unique.features = TRUE)
#write.table(HC,"Vis15_HN/C1/genes_mat.txt", sep="\t", quote = F, row.names = T)
geneMatHC<-"~/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Vis15_HN/C1/genes_mat.txt"
geneMatHC=fread(geneMatHC,sep="\t",header = T,stringsAsFactors = F,showProgress = T)
geneMatHC <- as.data.frame(geneMatHC) # force to data frame
rownames(geneMatHC) <- geneMatHC$GENE # make the name of rows GENEs
geneMatHC <- geneMatHC[,-1] # take ou
genesHC <- CreateSeuratObject(counts = HC, project = "HeadNeck_C")

colnames(geneMatHC)<-gsub("-1",".1",colnames(geneMatHC))
sample_combined_matPHNC<-rbind(geneMatHC,Sample_HMM_HC[,colnames(geneMatHC)])
sample_combined_matPHNC<-CreateSeuratObject(counts = sample_combined_matPHNC)
#sample_combined_matPHNCn<-subset(x = sample_combinedSB)

sample_combined_matPHNCn <- NormalizeData(sample_combined_matPHNC)
sample_combined_matPHNCn <- ScaleData(sample_combined_matPHNCn, features = rownames(genesHC))
#top2000<-head(VariableFeatures(genesPHNCn), 2000)
sample_combined_matPHNCn <- RunPCA(object = sample_combined_matPHNCn, reduction.name="pca_genes",features = rownames(genesHC))
sample_combined_matPHNCn <- FindNeighbors(object=sample_combined_matPHNCn, reduction="pca_genes", graph.name="gene_snn",dims=1:10, )
sample_combined_matPHNCn <- FindClusters(object=sample_combined_matPHNCn,graph.name="gene_snn",resolution=res_param)
sample_combined_matPHNCn <- RunUMAP(sample_combined_matPHNCn,reduction.name="umap_genes",reduction="pca_genes",dims=1:10,check_duplicates = F)

m<-read.csv("~/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Vis15_HN/C1/copykat/3subclones_BCs.txt", sep="\t", header=TRUE)
#m$BC<-gsub(".1","-1",m$BC)
sample_combined_matPHNCn<-sample_combined_matPHNCn[,m$BC]
rownames(m)<-m$BC
m <- m[match(colnames(sample_combined_matPHNCn), m$BC), ]
Idents(sample_combined_matPHNCn)<-m$clone
sample_combined_matPHNCn@meta.data$clone<-Idents(sample_combined_matPHNCn)
uTARs<-rownames(uTAR_seurat_HC)

sample_combined_matPHNCn<-SCTransform(sample_combined_matPHNCn, do.scale = FALSE,return.only.var.genes = FALSE)


# Check layers in RNA assay
print(names(sample_combined_matPHNCn@assays$RNA@layers))

# Get feature names from assay metadata or data layer
feature_names <- rownames(sample_combined_matPHNCn@assays$RNA@features)
if (is.null(feature_names) && "data" %in% names(sample_combined_matPHNCn@assays$RNA@layers)) {
    feature_names <- rownames(sample_combined_matPHNCn@assays$RNA@layers$data)
}
if (is.null(feature_names)) {
    stop("Cannot find feature names in RNA assay.")
}

# Assign row names to counts layer
if ("counts" %in% names(sample_combined_matPHNCn@assays$RNA@layers) && length(feature_names) == nrow(sample_combined_matPHNCn@assays$RNA@layers$counts)) {
    rownames(sample_combined_matPHNCn@assays$RNA@layers$counts) <- feature_names
    cat("Assigned row names to RNA assay counts layer.\n")
} else {
    stop("Feature names length does not match counts layer rows.")
}

# Verify uTARs
valid_uTARs <- uTARs[uTARs %in% feature_names]
cat("Number of valid uTARs:", length(valid_uTARs), "\n")

# Check non-zero expression
expr_matrix <- sample_combined_matPHNCn@assays$RNA@layers$counts[valid_uTARs, ]
feature_sums <- Matrix::rowSums(expr_matrix)
valid_uTARs <- valid_uTARs[feature_sums > 0]
cat("Number of uTARs with non-zero expression:", length(valid_uTARs), "\n")

# Set cluster identities
if ("pb_cluster" %in% colnames(sample_combined_matPHNCn@meta.data)) {
    Idents(sample_combined_matPHNCn) <- sample_combined_matPHNCn@meta.data$pb_cluster
}
table(Idents(sample_combined_matPHNCn))

# Run FindAllMarkers
Aall.markers <- FindAllMarkers(sample_combined_matPHNCn,
                               assay = "RNA",
                               slot = "data",
                               only.pos = TRUE,
                               min.pct = 0.05,
                               logfc.threshold = 0.05,
                               features = valid_uTARs,
                               test.use = "wilcox")
print(head(Aall.markers))
Atop20 <- Aall.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(sample_combined_matPHNCn, features = Atop20$gene, assay = "SCT") 


sample_filtered <- subset(sample_combined_matPHNCn, clone !="not.defined")
DoHeatmap(sample_filtered, features = Atop20$gene, assay = "SCT") 
