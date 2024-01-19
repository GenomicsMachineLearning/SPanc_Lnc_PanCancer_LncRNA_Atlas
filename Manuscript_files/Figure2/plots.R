## Phast Cons and MFE scatter plot
d<-read.csv("/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/Phastcons_MFE/HNC_phastcons_MFE_data.txt",sep=" ",header=FALSE)
highlight_names <- c("HOTAIR", "MALAT", "GAPDH", "RPL5", "HBA1","cuTAR213508","cuTAR213507","cuTAR100896")
ggplot(d, aes(x = V2, y = V3)) +
    geom_point(aes(color = ifelse(V1 %in% highlight_names, "red", "#b39eb5")), size = 1) +
    geom_text(data = subset(d, V1 %in% highlight_names), aes(label = V1), nudge_x = 0.1, nudge_y = 0.1) +
    scale_color_identity(guide = "none") +
    theme_minimal() + xlab("Conservation score") + ylab("Minimum Free Energy") + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

ggplot(d,aes(x=V2,y=V3))+geom_point()


#CPAT Coding potential BOXPLOT
d<-read.csv("/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/Phastcons_MFE/cpat_qqman.txt",sep="\t",header=TRUE)


