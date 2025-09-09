# library --------------
rm(list=ls());options(stringsAsFactors=FALSE)
library(tidyverse)
library(Seurat)
library(schard)
plan("multicore", workers = 24)
options(future.globals.maxSize = 250 * 1024^3)

# set path ------------------
output_dir <- '/scratch/project/stseq/Prakrithi/lnc_revision/xenium/mel/'
clean_dir <- file.path(output_dir,'clean-data')
input_dir <- '/QRISdata/Q4386/Xenium_20250319__042556__20250319_lncRNA_IO_Cancer/'
setwd(output_dir)
set.seed(111)

run_label_transfer <- function(sample, reference, assay="Spatial", ref_column="cell_type") {
    countsData <- sample@assays[[assay]]$counts
    Query_data <- SCTransform(sample, ncells = 2000, verbose = FALSE, vst.flavor = "v2", assay=assay)
    Query_data <- RunPCA(Query_data, verbose = FALSE, npcs=30)
    Query_data <- FindNeighbors(Query_data, reduction = "pca", dims = 1:30)
    Query_data <- FindClusters(Query_data, verbose = FALSE, resolution = 1.0)
    Query_data <- RunUMAP(Query_data, reduction = "pca", dims = 1:30)
    anchors <- FindTransferAnchors(
        reference = reference,
        query = Query_data,
        dims = 1:30,
        reference.reduction = "pca",
        normalization.method = "SCT",
        recompute.residuals = FALSE #, features=intersect(rownames(reference),rownames(Query_data))
    )
    predictions <- TransferData(
        anchorset = anchors,
        refdata = reference@meta.data[, ref_column],
        dims = 1:30
    )
    Query_data <- AddMetaData(Query_data, metadata = predictions)
    return(Query_data)
}


LoadXenium <- function (data.dir, fov = "fov", assay = "Xenium") {
    data <- ReadXenium(data.dir = data.dir, type = c("centroids", "segmentations"))
    segmentations.data <- list(centroids = CreateCentroids(data$centroids),
                               segmentation = CreateSegmentation(data$segmentations))
    coords <- CreateFOV(coords = segmentations.data, type = c("segmentation", "centroids"), molecules = data$microns, assay = assay)
    xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
    xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
    xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
    xenium.obj[[fov]] <- coords
    return(xenium.obj)
}


# reference with subset ----------------------
xenium_file <- 'output-XETG00114__0016987__Glioblastoma__20250319__042624'  # Specify your sample file name directly
sample_id = str_split_fixed(xenium_file, '__', n=4)[3]
xenium.obj <- LoadXenium(file.path(input_dir, xenium_file))
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 0)
samples <- xenium.obj

saveRDS(samples, file = file.path(clean_dir, 'xenium_GBM_loaded.rds'))

#samples<-readRDS("/QRISdata/Q4386/Xenium_20250319__042556__20250319_lncRNA_IO_Cancer/output-XETG00114__0016987__Colorectal__20250319__042625/CRC_loaded.rds")

# reference ----------------------
reference<-readRDS("/scratch/project_mnt/S0010/Prakrithi/lnc_revision/xenium/sc_ref/SCSP0000480_processed_downsampled.rds")
#reference <- readRDS("/scratch/project_mnt/S0010/Prakrithi/Vicki_CRC/SCP/CRC_sc_ref_GSE178341_crc10x_full_c295v4_submit.rds")
zero_count_cells <- colSums(GetAssayData(reference, assay='RNA', layer='counts')) == 0
reference <- reference[, !zero_count_cells]
reference <- SCTransform(reference, verbose = FALSE)
reference <- RunPCA(reference, verbose = FALSE)
reference <- FindNeighbors(reference, dims = 1:30, verbose = FALSE)
reference <- FindClusters(reference, verbose = FALSE)
reference <- RunUMAP(reference, dims = 1:30, verbose = FALSE)
reference2 <- reference

pdf("/scratch/project/stseq/Prakrithi/lnc_revision/xenium/mel/clean-data/GBM_ref.pdf", width = 10, height = 8)
DimPlot(reference2, group.by = "cell_type")
dev.off()

# label transform ----------------------
samples <- run_label_transfer(samples, reference = reference, ref_column = "cell_type", assay = "Xenium")
predicted_column <- str_detect(colnames(samples@meta.data), "^predict")
prediction <- samples@meta.data[, predicted_column]
write.csv(prediction, file.path(clean_dir, paste0(sample_id, "_LT_subsetted_ref.csv")))

saveRDS(samples, file = file.path(clean_dir, 'xenium_spatial_GBM_LT.rds'))

# save metadata --------------
write.csv(samples@meta.data, file = file.path(clean_dir, 'xenium_spatial_GBM_LT_meta.csv'), row.names = FALSE)
