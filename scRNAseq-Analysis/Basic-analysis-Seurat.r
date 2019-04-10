setwd("../Desktop")
library(Seurat)
library(dplyr)

#Formatting the data for Seurat
raw_data = read.table("geneExpDay0Day7.txt",header=TRUE)
raw_matrix = as.matrix(raw_data[,c(3:98)])
rownames(raw_matrix) = raw_data$SYMBOL
colnames(raw_matrix) = colnames(raw_data)[c(3:98)]

#Convert the data into a Seurat object
seurat_data = CreateSeuratObject(raw.data=raw_matrix)

#Normalization and centering the data
seurat_data = NormalizeData(object=seurat_data, normalization.method = "LogNormalize",scale.factor = 10000)
seurat_data = ScaleData(object=seurat_data)

#Run PC, only the first two PCs are computed in this case
seurat_data = RunPCA(object=seurat_data, pc.genes=raw_data$SYMBOL, do.print=TRUE, pcs.print=1:2, genes.print=5)

#Export the PCs in a table to generate a plot, this is feasible in Seurat with the PCAPlot function but it lacks graphical options.
PCA_data = FetchData(seurat_data,c("ident","PC1","PC2"))
col = vector(length=96)
col[PCA_data$ident=="ASP14.SC.day0"] = "red"
col[PCA_data$ident=="ASP14.SC.day7"] = "blue"
png("PCA.png", width=500, height=500)
plot(PCA_data$PC1,PCA_data$PC2, pch=16, col=col, ylab="PC2", xlab="PC1",main="PCA",ylim=c(-20,100))
legend("topleft", pch=16, col=c("red","blue"),legend=c("ASP14 day 0","ASP14 day 7"))
dev.off()

#Identification of the mislabelled cell, in this case ASP14.SC.day7_82 and modification of its status
which(PCA_data$PC1[PCA_data$ident=="ASP14.SC.day7"]<=-20)
PCA_data[PCA_data$ident=="ASP14.SC.day7",][34,]
ident = seurat_data@ident
ident[82] = "ASP14.SC.day0"
seurat_data = AddMetaData(seurat_data,metadata=ident,col.name="new.ident")
seurat_data = SetAllIdent(seurat_data, id='new.ident')

#Identifying genes characterizing the change from Day 0 to Day 7 with DE analysis, pvalue adjustment with Bonferroni
#The tobit test requires the VGAM package
Day7_Markers = FindMarkers(object=seurat_data, ident.1='ASP14.SC.day7', ident.2='ASP14.SC.day0',test.use='tobit')

#Set cutoff pvalue at 0.05
trimed_table = Day7_Markers[which(Day7_Markers$p_val_adj<0.05),]

#Generate a ranked gene list for GSEA, ranked on Pvalue or fold change
Pval_ranked_table = cbind(rownames(trimed_table), trimed_table[,5])
colnames(Pval_ranked_table) = c("Gene","Pval")
FC_trimed_table = trimed_table[order(abs(trimed_table$avg_logFC),decreasing=TRUE),]
FC_ranked_table = cbind(rownames(FC_trimed_table), FC_trimed_table[,2])
#If absolute fold change is preferred to fold change, use abs(FC_trimed_table[,2]) instead
colnames(FC_ranked_table) = c("Gene","FC")
write.table(Pval_ranked_table, "GSEA_table_2.rnk",col.names=TRUE,row.names=FALSE,quote=FALSE,sep='\t')
