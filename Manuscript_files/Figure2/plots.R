## Phast Cons and MFE scatter plot
d<-read.csv("/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/Phastcons_MFE/HNC_phastcons_MFE_data.txt",sep=" ",header=FALSE)
library(ggplot2)

# Assuming 'd' is your dataframe
# Create a new column 'color_group' to categorize the names into different colors
d$color_group <- with(d, ifelse(V1 %in% c("HOTAIR", "MALAT"), "orange",
                         ifelse(V1 %in% c("GAPDH", "RPL5", "HBA1"), "green",
                         ifelse(V1 %in% c("cuTAR213508", "cuTAR213507", "cuTAR100897", "cuTAR234975", "cuTAR170206"), "red", "purple"))))

# Plot
ggplot() +
  # Points for the rest (purple)
  geom_point(data = subset(d, color_group == "purple"), aes(x = V2, y = V3), color = "#b39eb5", size = 1) +
  # Points for highlighted names
  geom_point(data = subset(d, color_group != "purple"), aes(x = V2, y = V3, color = color_group), size = 1) +
  # Text labels for highlighted points, with a neutral color (e.g., black)
  geom_text(data = subset(d, V1 %in% highlight_names), 
            aes(x = V2, y = V3, label = V1, 
                vjust = ifelse(V1 %in% highlight_names, 1.5, 0), 
                hjust = ifelse(V1 %in% highlight_names, 0.5, 0.5)), 
            size = 3, color = "black") +  # Neutral color for text labels
  # Set colors manually with a colorblind-friendly green
  scale_color_manual(values = c("orange" = "orange", "green" = "#66CC00", "red" = "red", "purple" = "#b39eb5")) +
  theme_minimal() + 
  xlab("Conservation score") + 
  ylab("Minimum Free Energy") + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank()
  )



#CPAT Coding potential BOXPLOT
d<-read.csv("/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/Phastcons_MFE/cpat_qqman.txt",sep="\t",header=TRUE)


