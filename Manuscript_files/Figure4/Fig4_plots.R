################################# Moran's Alluvial ##############################################################################################################################
setwd("/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/morans/cancer_genes_moran/")
library(tidyverse)
library(igraph)
library(ggalluvial)
data <- read.csv("HNC_HH_3conf.txt",sep="\t")  

# Subset the data for counts greater than 50
data2 <- subset(data, count > 50)


# Calculate the total count for each hallmark
hallmark_counts <- data2 %>%
  group_by(hallmark) %>%
  summarize(total_count = sum(count))

# Reorder the hallmark factor levels based on the total count
hallmark_order <- hallmark_counts %>%
  arrange(desc(total_count)) %>%
  top_n(8, total_count) %>%
  pull(hallmark)

lncRNA_counts <- data2 %>%
  group_by(lncRNA) %>%
  summarize(total_count = sum(count))

lncRNA_order <- lncRNA_counts %>%
  arrange(desc(total_count)) %>%
  top_n(3, total_count) %>%
  pull(lncRNA)

data2 <- data2 %>%
  filter(hallmark %in% hallmark_order)

data2$hallmark <- factor(data2$hallmark, levels = rev(hallmark_order))
# Reorder the lncRNA factor levels based on the count
#data2$lncRNA <- reorder(data2$lncRNA, -data2$count)
data2$lncRNA <- factor(data2$lncRNA, levels = rev(lncRNA_order))

#color_palette<-c("#2166AC","#B2182B","#D6604D","#F4A582","#FDDBC7","#e8ddb7","#999999","#D1E5F0","#92C5DE","#4393C3","#053061","#67001F","brown")
color_palette<-c("#2166AC","#B2182B","#D6604D","#F4A582","#FDDBC7","#e8ddb7","#999999","#D1E5F0")

#color_paletteh<-c("#2166AC","#B2182B","#D6604D","#F4A582","#FDDBC7","#e8ddb7","#999999","#D1E5F0") #for 8 hallmarks
# Create the alluvial plot
ggplot(data2, aes(y = count, axis1 = lncRNA, axis2 = hallmark)) +
  geom_alluvium(aes(fill = hallmark, alpha = count),curve_type = "cubic", width = 2/4, knot.pos = 1, reverse = FALSE) +
  guides(fill = "none") +
  geom_stratum(alpha = 0, width = 2/4, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum),  fontface = "bold"), reverse = FALSE, size = 5, hjust = 0.5) +
  scale_x_continuous(breaks = 1:2, labels = c("lncRNA", "Hallmark")) +
  scale_fill_manual(values = rev(color_palette))+
  theme(axis.text = element_text(size = 10)) #10 for top uTARs and width = 2/4
# save pdf landscape 9x8

