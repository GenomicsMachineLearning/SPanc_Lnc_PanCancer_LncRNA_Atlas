# Load libraries
lapply(c("dplyr","igraph","Seurat","ggplot2","ggrepel","enrichR","data.table","caTools","scales","stringr","cluster"), library, character.only = TRUE)

#library(rolypoly)
#set parameters
nDiffFeat<-5
cells<-2
feats<-100
nFeatMin<-200
nFeatMax<-2500
pcaDims<-30
res_param<-0.2

setwd("~/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Vis15_HN/C1/")
#setwd("~/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Alex_breast/s4")

#uTARs
noDirMat38_PHNC<-("uTARs_uniq_mat_gc43.txt")
noDirMat38_PHNC<-("HNC_renamed_final_mat.txt")
#noDirMat38_PHNC<-("../vis/C/uniq_uTARs_mat")
Sample_HMM_PHNC=fread(noDirMat38_PHNC,sep="\t",header = TRUE,stringsAsFactors = FALSE,showProgress = TRUE)
rownames(Sample_HMM_PHNC)<-tolower(rownames(Sample_HMM_PHNC))
Sample_HMM_PHNC <- as.data.frame(Sample_HMM_PHNC) # force to data frame
rownames(Sample_HMM_PHNC) <- Sample_HMM_PHNC$GENE # make the name of rows GENEs
Sample_HMM_PHNC <- Sample_HMM_PHNC[,-1] # take out first column
uTAR_seurat_PHNC<-CreateSeuratObject(counts = Sample_HMM_PHNC, project = "HNC",min.cells = 3, assay = "Spatial")


#GENES
testPHNC <-  Read10X(data.dir ="filtered_feature_bc_matrix/", gene.column = 2, unique.features = TRUE)
geneMatPHNC<-"genes_mat.txt"
#testPHNC <-  Read10X(data.dir ="../vis/C/outs/filtered_feature_bc_matrix/", gene.column = 2, unique.features = TRUE)
#geneMatPHNC<-"../vis/C/genes_mat.txt"
geneMatPHNC=fread(geneMatPHNC,sep="\t",header = TRUE,stringsAsFactors = FALSE,showProgress = TRUE)
geneMatPHNC <- as.data.frame(geneMatPHNC) # force to data frame
rownames(geneMatPHNC) <- geneMatPHNC$GENE # make the name of rows GENEs
geneMatPHNC <- geneMatPHNC[,-1] # take ou
genesPHNC <- CreateSeuratObject(counts = testPHNC, project = "HNC", assay="Spatial")

genesPHNCn <- NormalizeData(genesPHNC)
genesPHNCn <- ScaleData(genesPHNCn, features = rownames(genesPHNCn))
genesPHNCn <- FindVariableFeatures(genesPHNCn, selection.method = "vst", nfeatures = 2000)

# combine genes and uTARs
sample_combined_matPHNC<-rbind(geneMatPHNC,Sample_HMM_PHNC[,colnames(geneMatPHNC)])
sample_combined_matPHNC<-CreateSeuratObject(counts = sample_combined_matPHNC, assay = "Spatial")

sample_combined_matPHNCn <- NormalizeData(sample_combined_matPHNC)
sample_combined_matPHNCn <- ScaleData(sample_combined_matPHNCn, features = rownames(sample_combined_matPHNCn))
#top2000<-head(VariableFeatures(genesPHNCn), 2000)
sample_combined_matPHNCn <- RunPCA(object = sample_combined_matPHNCn, reduction.name="pca_genes",features = VariableFeatures(genesPHNCn))
sample_combined_matPHNCn <- FindNeighbors(object=sample_combined_matPHNCn, reduction="pca_genes", graph.name="gene_snn",dims=1:10, )
sample_combined_matPHNCn <- FindClusters(object=sample_combined_matPHNCn,graph.name="gene_snn",resolution=res_param)
sample_combined_matPHNCn <- RunUMAP(sample_combined_matPHNCn,reduction.name="umap_genes",reduction="pca_genes",dims=1:10,check_duplicates = F)

# Rename clusters based on oufrom scType
Idents(sample_combined_matPHNCn)<-Idents(spatial_data)
DimPlot(sample_combined_matPHNCn, reduction = "umap_genes")


######## Load libraries for WGCNA #####
# single-cell analysis package
library(Seurat)
# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA) # co-expression network analysis packages:
#devtools::install_github('smorabit/hdWGCNA', ref='dev')
library(hdWGCNA)
theme_set(theme_cowplot()) # using the cowplot theme for ggplot
set.seed(12345) # set random seed for reproducibility

#set up visium object
c<-read.csv("spatial_anndata",sep="\t", header=TRUE) #read spatial data tissue,row,col,imagerow,imagecol
#c$BC<-gsub('.1','-1',c$BC)
colnames(c)<-c("X","tissue","row","col","imagerow","imagecol")
rownames(c) = c$X
c<-c[,-1]
sample_combined_matPHNCn@meta.data<-c #add to metadata
sample_combined_matPHNCn@meta.data$cell_type<-Idents(sample_combined_matPHNCn)

##rename clusters based on Morans
m<-read.csv("/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/morans/cancer_genes_moran/forWGCNA_labels.txt", sep="\t", header=TRUE)
Idents(sample_combined_matPHNCn)<-m$GAPDH.chr18_54606899_54607299_._2526_0
sample_combined_matPHNCn$Morans<-Idents(sample_combined_matPHNCn)


# read image data from image file
img<-Read10X_Image(
  "spatial/",
  image.name = "tissue_lowres_image.png",
  filter.matrix = TRUE,
)
sample_combined_matPHNCn@images$HNC =  img # add image data to slot
sample_combined_matPHNCn@images$HNC@key <- "HNC_" #key SHOULD be same as the image slot name
sample_combined_matPHNCn@images$HNC@assay <- "Spatial"
rownames(sample_combined_matPHNCn@images$HNC@coordinates)<-gsub("-1",".1",rownames(sample_combined_matPHNCn@images$HNC@coordinates))

saveRDS(sample_combined_matPHNCn,"sample_combined_matPHNCn.rds")


####setting up for WGCNA
HNC_wgcna<- SetupForWGCNA(sample_combined_matPHNCn, gene_select = "features", # the gene selection approach
                          features = rownames(sample_combined_matPHNCn), # fraction of cells that a gene needs to be expressed in order to be included
                          wgcna_name = "HNC" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seu_wgcna <- MetaspotsByGroups(
  HNC_wgcna,
  group.by = c("cell_type"),
  ident.group = "cell_type",
  assay = 'Spatial',
  slot = 'counts'
)

#### construct metacells  in each group based on morans
seu_wgcna <- MetaspotsByGroups(
  HNC_wgcna,
  group.by = c("Morans"),
  ident.group = "Morans",
  assay = 'Spatial',
  slot = 'counts'
)
# normalize metacell expression matrix:
seu_wgcna <- NormalizeMetacells(seu_wgcna)
m_obj <- GetMetacellObject(seu_wgcna) 
m_obj

## Co-expression network analysis

## to include all spots
#seu_wgcna_all  <- SetDatExpr(seu_wgcna,group.by=NULL,group_name = NULL) 

# set-up for cancer cells
seu_wgcna <- SetDatExpr(
  seu_wgcna,
  group_name = "Cancer cells", # the name of the group of interest in the group.by column
  group.by='cell_type', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'Spatial', # using RNA assay
  slot = 'data' # using normalized data
)


# set-up for Morans HH
seu_wgcna <- SetDatExpr(
  seu_wgcna,
  group_name = "HH", # the name of the group of interest in the group.by column
  group.by='Morans', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'Spatial', # using RNA assay
  slot = 'data' # using normalized data
)

# select soft-power threshold
# Test different soft powers:
seu_wgcna <- TestSoftPowers(
  seu_wgcna,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seu_wgcna)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seu_wgcna)
head(power_table)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seu_wgcna, soft_power=9, #set soft_power based on plot # 8 for seu_wgcna_all # 12 for HH Morans
  setDatExpr=FALSE,
  tom_name = 'Cancer cells' # name of the topoligical overlap matrix written to disk
)  #Takes long

seurat_obj <- ConstructNetwork(
  seu_wgcna, soft_power=14, #set soft_power based on plot # 14 for HH chr21 SAT1 # 12 for HH Morans FGFR3
  setDatExpr=FALSE,
  tom_name = 'HH' # name of the topoligical overlap matrix written to disk
) 

PlotDendrogram(seurat_obj, main='Tumor HNC hdWGCNA Dendrogram')
all_dend<-PlotDendrogram(seurat_obj_all, main='All HNC hdWGCNA Dendrogram')
PlotDendrogram(seurat_obj, main='Morans HH hdWGCNA Dendrogram')



## Inspecting the TOM
TOM <- GetTOM(seurat_obj)

## COMPUTE hormonized module eigengenes - only when integrated seurat objects are used
# need to run ScaleData first or else harmony throws an error:
#seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars=NULL
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)


# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cell_type', group_name = 'Cancer cells', sparse=FALSE
)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'Morans', group_name = 'HH', sparse=FALSE
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "HNC-M"
)

# plot genes ranked by kME for each module
pkME <- PlotKMEs(seurat_obj, ncol=5)
pkME

##module network plots
ModuleNetworkPlot(seurat_obj)

# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'moi' #'all
) # gets Killed

#umap network
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)

#
# get the module assignment table:
modules <- GetModules(seurat_obj)
grep("cuTAR213507",rownames(modules)) # grep uTAR of interest cuTAR170206
modules_all <- GetModules(seurat_obj_all)
write.table(modules,"modulesn.txt", sep="\t", row.names = TRUE, quote = FALSE)

# show the first 6 columns:
head(modules[,1:6])

#### Differential Module Expression
group1 <- seurat_obj@meta.data %>% subset(cell_type == 'Cancer cells') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(cell_type == 'Naive CD4+ T cells') %>% rownames


group1 <- seurat_obj@meta.data %>% subset(Morans == 'HH') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(Morans = c('LL','LH','HL','Non-Significant')) %>% rownames
group2 <- group2[!(group2 %in% group1)]

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='HNC'
)

PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  wgcna_name = 'HNC'
)


#### MODULE ENRICHMENT ANALYSIS ###
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 500 # number of genes per module to test
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)


top20_mods<-DMEs[1:20,]
top20_modules<-top20_mods$module
uTAR_mods<- c("Tumor-M9","Tumor-M49","Tumor-M20")

# enrichr dotplot
EnrichrDotPlot_custom(
  seurat_obj,
  mods = "HNC-M14",  #M14 for chr21# use all modules "all" (this is the default behavior) or specify module name
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=10 # number of terms for each module
)


## corr
corr<-ModuleCorrelogram(seurat_obj, features="MEs")
write.table(corr$corr,"corr.txt",sep="\t", quote = FALSE, row.names = TRUE)



all.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top50 <- all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
DoMultiBarHeatmap(seurat_obj, assay = 'RNA', features = top_10$gene, group.by='Idents', additional.group.by = 'orig.ident')

alluTAR.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1, features = rownames(uTAR_seurat_PHNC))
top50_uTAR <- alluTAR.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
write.table(top50_uTAR,"HNC_DE_uTARs.txt",sep="\t",row.names = FALSE, quote=FALSE)


# plot using ggplot2
ggplot(df, aes(x = Position, y = LogFC, color = CellType)) + 
  geom_point() + 
  labs(x = "Gene Position", y = "LogFC", color = "Cell Type")

ggplot(d, aes(x = gene, y = avg_log2FC, color = cluster)) + 
  geom_point() + 
  labs(x = "Gene ", y = "LogFC", color = "Cell Type")

# plot using ggplot2
ggplot(d, aes(x = gene, y = avg_log2FC, fill = feature)) + 
  geom_point() + 
  facet_wrap(~ cluster, ncol = 1) +
  labs(x = "Gene ", y = "LogFC", fill = "Feature")



modnetgenes<-c("cuTAR170206","ZFP36","NDRG1","ALDOA","PFKFB3","IFNGR2","UPP1","NAMPT","NFKBIA","ADM","VEGFA","ERRFI1","IGFBP3","DDIT4","MKNK2","IL1A","SLC16A3","HS3ST1","NUPR1")
cur <- cur[, c("gene_name", "kME_HNC_M14")]

# Subset 'cur' based on matching genes in 'modnetgenes'
cur <- subset(modules, module == "HNC-M14")
cur <- cur[, c("gene_name", "kME_HNC-M14")]
# Subset 'cur' based on matching genes in 'modnetgenes'
cur <- cur[cur$gene_name %in% modnetgenes, ]

