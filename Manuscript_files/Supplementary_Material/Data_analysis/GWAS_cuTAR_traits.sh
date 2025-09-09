# On terminal
bedtools intersect -a /QRISdata/Q4386/all_cuTARs/all_cuTARs_refined_after_gencode47_clean_v5.bed -b ../GWAS_catalog_All_associations_v1.02.bed -wo | awk -F"\t" '{print $4"\t"$10}' | sed 's/ /_/g' > utar_trait.txt
awk 'BEGIN { FS = OFS = "\t" } { split($2, arr, /,/); $2 = arr[1] } 1' utar_trait.txt > utar_trait2.txt

terms=( $(cut -f 2 utar_trait2.txt | sort -u) )   # Store the unique terms from column 2 in an array
for term in "${terms[@]}"; do
  count=$(grep -w "$term" utar_trait2.txt | cut -f 1 | sort | uniq | wc -l)
  echo "$term $count"
done > uTAR_trait_counts.txt

sort -k2,2nr uTAR_trait_counts.txt | head 
awk '/cancer|carcinoma|rcoma|noma|toma/' utar_trait2.txt | awk  '{print $1}' | sort | uniq | wc -l --> GWAS Cancer cuTARs
1144
awk '/cancer|carcinoma|rcoma|noma|toma/' uTAR_trait_counts.txt > uTAR_cancer_trait_counts.txt
# R plot

## in R
# Assuming your data frame `d` has columns: term and count
library(ggplot2)
d<-read.csv("uTAR_cancer_trait_counts.txt",sep=" ", header=FALSE)
# Reorder terms based on count descending
d$V1 <- factor(d$V1, levels = d$V1[order(d$V2, decreasing = FALSE)])

# Create the barplot
ggplot(d, aes(x = V1, y = V2)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip for better readability if many terms
  labs(x = "Term", y = "Count", title = "Term Counts (Descending)") +
  theme_minimal()
