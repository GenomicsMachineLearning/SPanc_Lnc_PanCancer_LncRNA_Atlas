
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
