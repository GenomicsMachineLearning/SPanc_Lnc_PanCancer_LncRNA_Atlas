library(ggpubr)
my_colors <- c("Pre" = "#967bb6", "Post" = "#fa8072")
p<-read.csv("/Users/ritanya/Downloads/Acral/A3_Treatment/pseudo/new/Random10_pseudo_bulk.txt",sep="\t",header=TRUE)
p$Condition <- factor(p$Condition, levels = c("Pre", "Post"))
p <- subset(p, Cell_type != "Unknown")
uTAR8337<-p[p$gene %in% "uTAR8337",] # good
uTAR2929<-p[p$gene %in% "uTAR2929",] # good

p1<-ggplot(uTAR8337, aes(x = reorder(Cell_type, -Expression), y = Expression, fill = Condition)) +
  geom_boxplot() +
  scale_fill_manual(values = my_colors) +  # Set manual colors
  theme_minimal() +  # Minimal theme
  theme(
    plot.background = element_rect(fill = "white"),  # White background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove panel border
    axis.line = element_line(),  # Keep axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "bottom",  # Move legend to the bottom
    axis.ticks = element_line()  # Keep tick lines
  ) +
  labs(x = "Cell Type", y = "Expression") +  # Customize axis labels
  guides(fill = guide_legend(reverse = TRUE)) +  # Reverse legend order
  stat_compare_means(aes(group = Condition))

#combined facet
combined_df <- rbind(uTAR2929,uTAR8337)
combined_df$Source <- factor(rep(c("uTAR8337", "uTAR2929"), each = nrow(uTAR8337)))
ggplot(combined_df, aes(x = reorder(Cell_type, -Expression), y = Expression, fill = Condition)) +
  geom_boxplot() +
  scale_fill_manual(values = my_colors) +  # Set manual colors
  theme_minimal() +  # Minimal theme
  theme(
    plot.background = element_rect(fill = "white"),  # White background
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_blank(),  # Remove panel border
    axis.line = element_line(),  # Keep axis lines
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "bottom",  # Move legend to the bottom
    axis.ticks = element_line()  # Keep tick lines
  ) +
  labs(x = "Cell Type", y = "Expression") +  # Customize axis labels
  guides(fill = guide_legend(reverse = TRUE)) +  # Reverse legend order
  facet_wrap(~ Source, ncol = 1, strip.position = "right")
