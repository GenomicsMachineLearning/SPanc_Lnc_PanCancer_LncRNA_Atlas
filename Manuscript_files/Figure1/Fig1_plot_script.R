########################  Fig1 barplot ############################
data<-read.csv("lncAtlas_barplot_data", sep="\t", header=TRUE)
d<-pivot_longer(
  data = data,
  cols = -Sample,
  names_to = "LncRNAs",
  values_to = "No. of lncRNAs"
)
d$Sample <- factor(d$Sample,levels=unique(d$Sample))

ggplot(d, aes(x = Sample, y = d$`No. of lncRNAs`, fill = LncRNAs)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_grey() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom"
  ) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + ylab("No. of LncRNAs")

ggsave("stats_barplot_600dpi.png", plot = p, width = 8.3, height = 5, units = "in", dpi = 600)

# or
#hex codes google colors:  https://spreadsheet.dev/how-to-get-the-hexadecimal-codes-of-colors-in-google-sheets
custom_colors <- c("Novel.uTARs" = "#46bdc6", "lncExpDB" = "#f9cb9c", "FANTOM" = "#b6d7a8", "FANTOM...lncExpDB" = "#ffe599")
p<-ggplot(d, aes(x = Sample, y = `No. of lncRNAs`, fill = LncRNAs)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = custom_colors) +  # Set custom colors here
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom"
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
  ylab("No. of LncRNAs")


################################# Fig 1 Relative gene- utar counts violin ######################
library(viridis)
v<-read.csv("/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/relative_featurecounts_violin_new", sep="\t", header=TRUE)
pastel_colors <- c(
  "#FFD700", "#FFA07A", "#90EE90", "#87CEFA", "#FFC0CB",
           "#FFDAB9", "#98FB98", "#DDA0DD", "#ADD8E6", "#FFB6C1"
)
purples <- c("#F1E4F6", "#E3D0EE", "#D5BCE6", "#C7A8DE", "#B995D7", "#AA81CF", 
                         "#9C6DC7", "#8E59BF", "#8044B7", "#7220AF", "#631DA0", "#541B91", "#461883")
                         
v$Sample <- factor(v$Sample,levels=unique(v$Sample))
ggplot(v, aes(y=as.factor(Samples.1), x=as.numeric(log(Relative_counts)), fill=Samples.1)) + geom_violin() +theme_classic() + xlab("Relative Gene-uTAR counts (log)") + ylab("Samples")+scale_fill_manual(values = pastel_colors)
ggplot(v, aes(y = as.factor(Sample), x = as.numeric(log), fill = Sample)) +
  geom_violin() +
  theme_classic() +
  xlab("Relative Gene-uTAR counts (log)") +
  ylab("Samples") +
  scale_fill_manual(values = inferno(length(unique(v$Sample)))) +
  theme(
    axis.text = element_text(size = 12),   # Adjust the size of tick labels
    axis.title = element_text(size = 14)  # Adjust the size of axis labels
  )
