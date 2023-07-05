#merged_df is your raw gene expression matrix
#col -> samples; row -> genes
edgeR_meta<-combined@meta.data[,c("gene_cell_types","orig.ident")]
write.table(edgeR_meta,"edgeR_metadata.txt",sep="\t",quote=FALSE,row.names=TRUE)

#sed '1s/.1/.1_1/g' A1_gene_uTARs_final_mat.txt > edgeR/A1_mat.txt
#sed '1s/.1/.1_2/g' A2_gene_uTARs_final_mat.txt > edgeR/A2_mat.txt
#sed '1s/.1/.1_3/g' A5_gene_uTARs_final_mat.txt > edgeR/A5_mat.txt
#sed '1s/.1/.1_4/g' A6_gene_uTARs_final_mat.txt > edgeR/A6_mat.txt
#sed '1s/.1/.1_5/g' C1_gene_uTARs_final_mat.txt > edgeR/C1_mat.txt
A1<-read.csv("A1_mat.txt",sep="\t",header=TRUE)
A2<-read.csv("A2_mat.txt",sep="\t",header=TRUE)
A5<-read.csv("A5_mat.txt",sep="\t",header=TRUE)
A6<-read.csv("A6_mat.txt",sep="\t",header=TRUE)
C1<-read.csv("C1_mat.txt",sep="\t",header=TRUE)

#merged_df <- merge(A1, A2, A5, A6, C1, by = 'GENE', all = TRUE)
merged_df <- cbind(A1, A2[, -1], A5[, -1], A6[, -1], C1[, -1])
write.table(merged_df,"combined_mat_raw.txt",sep="\t",quote=FALSE, row.names=FALSE)

###
library(edgeR)
merged_df<-read.csv("combined_mat_raw.txt",sep="\t",header=TRUE) 
rownames(merged_df)<-merged_df$GENE
merged_df<-merged_df[,-1]
#Read metadata
meta<-read.csv("edgeR_metadata.txt",sep="\t",header=TRUE) 
is_mito<-grepl('mt-',row.names(merged_df))
is_ribo<-grepl('Rps|Rpl',row.names(merged_df))
keep_gene <- !is_mito & !is_ribo
#y <- DGEList(merged_df, samples=samples,group=group) # DE between types: YES AND NO
y <- DGEList(merged_df, samples=meta$orig.ident,group=meta$gene_cell_types) # DE between types: YES AND NO
#y$samples$group <- group
#keep <- filterByExpr(y, group = meta$gene_cell_types, min.count = 2, min.total.count = 10) #set your own threshold
#y <- y[keep,]
patients <- factor(meta$orig.ident)
y <- calcNormFactors(y) #takes time
design <- model.matrix(~y$samples$group+patients)
y <- estimateDisp(y, design, robust=TRUE) #takes time
fit <- glmQLFit(y, design, robust=TRUE) #probably takes time
res <- glmTreat(fit, coef=ncol(fit$design), lfc=1) #takes time
summary(decideTests(res))
res$table$regulate <- recode(as.character(decideTests(res)),"0"="Normal","1"="Up","-1"="Down")
de_group_edgeR <- res$table[order(res$table$PValue),]
table(decideTests(res))
res <- topTags(res,n = nrow(y))
res$table$regulate <- "Normal"
res$table$regulate[res$table$logFC>0 & res$table$FDR<0.05] <- "Up"
res$table$regulate[res$table$logFC<0 & res$table$FDR<0.05] <- "Down"
de_group_edgeR <- res$table[order(res$table$FDR),]
table(de_group_edgeR$regulate)
saveRDS(de_group_edgeR,"de_group_edgeR.rds")
write.table(table(de_group_edgeR$regulate),"de_edgeR_results.txt",sep="\t",quote=FALSE,row.names = TRUE)


