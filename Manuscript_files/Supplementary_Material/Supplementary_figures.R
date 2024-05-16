## Supplementary Fig. 4 Upset plots
ibrary(UpSetR)
library(ComplexHeatmap)
df<-read.csv("/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/pan_cancer_upset.txt", sep="\t", header=TRUE)
df<-read.csv("/Users/ritanya/Downloads/pan_cancer_upset_no_lung.txt", sep="\t", header=TRUE)
m = make_comb_mat(df, top_n_sets = 11)
mh<-m[comb_degree(m) >= 9]
ml<-m[comb_degree(m) ==1]

UpSet(m, top_annotation = upset_top_annotation(m, add_numbers = TRUE),
      right_annotation = upset_right_annotation(m, add_numbers = TRUE))

###############################################################################################################################
## Supplementary Figure 10a Downsampled nanopore UMI-BC plot
d<-read.csv("UMII_BC_sampling",sep="\t",header=TRUE)
df <- d %>% filter(No..of.UMIs < 1000)
ggplot(df, aes(x = as.factor(No..of.UMIs), y = No.of.uTARs)) + geom_boxplot()+  theme_minimal() +labs(x = "No. of UMIs per spot", y = "No.of uTARs")

###############################################################################################################################
## Supplementray Fig. 11
# Dot plot for % overlap of cuTAR expressing spots across Visium and long read
cp <- as.data.frame(read_excel("/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/nano_10x_overlap.xlsx", sheet = 8))
cp$Sample <- factor(cp$Sample, levels = c("HNC","SCC","BCC"))
ggplot(cp, aes(x = Sample, y = cuTAR, size = -log(p + 0.0000001), color = overlap_percent, label = paste("Overlap:", round(overlap_percent), "%"))) +
  geom_point() +
  geom_text(aes(label = paste("", round(overlap_percent),"%")), color = "black", size = 2.5) + 
  scale_size_continuous(range = c(1, 10), trans = "log") +
  labs(
    title = "%Overlap between Visium and ONT",
    x = "Sample",
    y = "cuTAR",
    size = "-log(P-Value)",
    color = "Percentage"
  ) +
  scale_color_gradient2(low = "yellow", mid = "orange", high = "red", midpoint = 15) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###############################################################################################################################
## Supplementray Fig. 17a
# Dot plot for % co-expressed spots across cell types
df<-read.csv("/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Alex_breast/s4/scType_clusters.txt", sep="\t", header=TRUE)
calculate_percentage_with_p_value <- function(df, column_name) {
  df %>%
    group_by(cell_type) %>%
    summarise(
      percentage_HH = mean(.data[[column_name]] == "HH", na.rm = TRUE),
      percentage_LL = mean(.data[[column_name]] == "LL", na.rm = TRUE),
      percentage_LH = mean(.data[[column_name]] == "LH", na.rm = TRUE),
      percentage_HL = mean(.data[[column_name]] == "HL", na.rm = TRUE),
      percentage_NS = mean(.data[[column_name]] == "Non-Significant", na.rm = TRUE),
      p_valueHH = ifelse(length(unique(.data[[column_name]])) > 1, binom.test(sum(.data[[column_name]] == "HH"), length(.data[[column_name]]))$p.value, 1),
      p_valueLL = ifelse(length(unique(.data[[column_name]])) > 1, binom.test(sum(.data[[column_name]] == "LL"), length(.data[[column_name]]))$p.value, 1),
      p_valueLH = ifelse(length(unique(.data[[column_name]])) > 1, binom.test(sum(.data[[column_name]] == "LH"), length(.data[[column_name]]))$p.value, 1),
      p_valueHL = ifelse(length(unique(.data[[column_name]])) > 1, binom.test(sum(.data[[column_name]] == "HL"), length(.data[[column_name]]))$p.value, 1),
      p_valueNS = ifelse(length(unique(.data[[column_name]])) > 1, binom.test(sum(.data[[column_name]] == "Non-Significant"), length(.data[[column_name]]))$p.value, 1)
    )
}
# Calculate the percentage with p-value for the "BIRC5" column
result_BIRC5 <- calculate_percentage_with_p_value(df, "BIRC5")
result_ELAVL1 <- calculate_percentage_with_p_value(df, "ELAVL1")
result_TARDBP <- calculate_percentage_with_p_value(df, "TARDBP")
result_BIRC5_ELAVL1 <- calculate_percentage_with_p_value(df, "BIRC5_ELAVL1")

# Print the result
print(result_BIRC5)
print(result_ELAVL1)
print(result_TARDBP)
print(result_BIRC5_ELAVL1)

# Calculate the percentage of "HH" for each cell type under the "BIRC5" column
result_df <- df %>%
  filter(BIRC5 == "HH") %>%
  group_by(cell_type) %>%
  summarise(
    percentage_HH = n() / nrow(.)
  )
# Print the result
print(result_df)

# Initialize an empty data frame to store results
result_df <- data.frame(cell_type = character(), value = character(), percentage = numeric(), p_value = numeric(), stringsAsFactors = FALSE)
# Loop through each value in BIRC5
for (value in unique_values) {
  # Loop through each cell type
  for (cell_type in unique_cell_types) {
    # Calculate percentage for the specific combination
    percentage <- sum(df$cell_type == cell_type & df$BIRC5 == value) / sum(df$BIRC5 == value)
    
    # Run Fisher's exact test for the specific combination
    contingency_table <- table(df$cell_type == cell_type, df$BIRC5 == value)
    fisher_test_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
    
    # Append results to the data frame
    result_df <- rbind(result_df, data.frame(cell_type = cell_type, value = value, percentage = percentage, p_value = fisher_test_result$p.value))
  }
}
#print(result_df)
result_df1 <- result_df %>%
  filter(value != "Non-Significant" & cell_type != "Unknown")

p1<-ggplot(result_df1, aes(x = cell_type, y = value, size = -log(p_value+0.0000001), color = percentage*100)) +
  geom_point() +
  scale_size_continuous(range = c(1, 10), trans = "log") +
  labs(
    title = "cuTAR215705 ~ BIRC5",
    x = "Cell Type",
    y = "Spatial Autocorrelation",
    size = "-log(P-Value)",
    color = "Percentage"
  ) +
  theme_minimal() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 35) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###############################################################################################################################


p1 + p2 + plot_layout(ncol = 2, heights = c(3, 3))

