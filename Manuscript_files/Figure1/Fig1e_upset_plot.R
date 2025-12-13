# Load VennDetail package
# BiocManager::install("VennDetail")
library(VennDetail)
options(repr.plot.height = 10, repr.plot.width = 15)

# Read existing sets
all       <- read.csv("/QRISdata/Q4386/all_cuTARs/cuTARs_final_list_v5_new", header = FALSE)
full      <- read.csv("characterization_re/split/correct/overlap_outputs_final/full_transcript_cuTARs_list", header = FALSE)
last      <- read.csv("characterization_re/split/correct/overlap_outputs_final/last_exon_cuTARs_list", header = FALSE)
internal  <- read.csv("characterization_re/split/correct/overlap_outputs_final/internal_exon_cuTARs_list", header = FALSE)
elncrna   <- read.csv("characterization_re/split/correct/elncRNA_1000bp_upstream_list", header = FALSE)
none      <- read.csv("characterization_re/split/correct/overlap_outputs/no_evidence_cuTARs", header = FALSE)
gencode47 <- read.csv("/scratch/project/stseq/Prakrithi/lnc_revision/In_gencode_cuTARs_list_new", header = FALSE)

# Read TSS and PAS evidence sets (newly added)
tss       <- read.csv("characterization_re/split/correct/overlap_outputs_final/cutars_with_TSS.txt", header = FALSE)
pas       <- read.csv("characterization_re/split/correct/overlap_outputs_final/cutars_with_PAS.txt", header = FALSE)

# Create VennDetail object with all sets including TSS and PAS
ven <- venndetail(list(
  All_cuTARs        = all$V1,
  Complete_evidence = full$V1,
  Last_Exon         = last$V1,
  Internal_Exon     = internal$V1,
  elncRNA           = elncrna$V1,
  TSS_evidence      = tss$V1,
  PAS_evidence      = pas$V1,
  No_evidence       = none$V1,
  Gencode47         = gencode47$V1
))

# Plot vennpie
plot(ven, type = "upset")
# Install if needed
# install.packages("ComplexUpset")

library(ComplexUpset)
library(dplyr)

# Read all your cuTAR sets
input_list <- list(
  All_cuTARs        = all$V1,
  Complete_evidence = full$V1,
  Last_Exon         = last$V1,
  Internal_Exon     = internal$V1,
  TSS_evidence      = tss$V1,
  PAS_evidence      = pas$V1,
  Gencode47         = gencode47$V1
)

# Make a presence/absence matrix
all_tars <- unique(unlist(input_list))
binary_matrix <- data.frame(cuTAR = all_tars)

for (set_name in names(input_list)) {
  binary_matrix[[set_name]] <- as.integer(binary_matrix$cuTAR %in% input_list[[set_name]])
}

# Plot using ComplexUpset with black bars sorted by size
p<-upset(
  binary_matrix,
  intersect = names(input_list),
  sort_intersections_by = "cardinality",   # <- This sorts black bars!
  sort_intersections = "descending",       # <- In descending order
  name = "cuTARs",
  base_annotations = list(
    'Intersection size' = intersection_size()
  )
)
p
