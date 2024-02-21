## MSCs RNA-seq data analysis. 
# Samples from tooths: PDLSCs, DPSCs, GFs. 
# Cell status: undifferentiated and 10 days after osteogenic differentiation.

# Code was created by Polina Kuchur
install.packages("readxl")
library(readxl)
library(org.Hs.eg.db)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(edgeR)
library(ggplot2)
library(plotly)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(ReactomePA)
library(clusterProfiler)
library(VennDiagram)


setwd("C:/Users/Юлия/Desktop/Анфиса/Уроки/ITMO/R_praktics/ligament_transcriptome")

count_matrix <- as.matrix(read.csv("./ligament_transcriptomes_total.csv", 
                                   sep = ",", row.names = 1))
head(count_matrix, 2)

#metadata (information about samples)
snames <- colnames(count_matrix)
status <- substr(snames, 18, nchar(snames) - 2)
status <- as.factor(status)

coldata <- data.frame(
  sample = c(colnames(count_matrix)),
  replicates = substr(snames, 1, nchar(snames) - 2),
  status = status,
  batch = substr(snames, 1, 16),
  row.names = "sample")

#check the order of samples in metadata and count_matrix
all(rownames(coldata) %in% colnames(count_matrix))


# prepare data for DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~batch + status)

dds <- collapseReplicates(dds, dds$replicates)

# set the reference condition (undif)
dds$status <- relevel(dds$status, ref = "48h.undif")

#filtration of low-counted genes
dds <- dds[rowSums(counts(dds)) >= 10,]
nrow(dds)  #считает кол-во строчек в аргументе

## Add gene symbols to gene properties
genes <- rownames(dds)
symbols <- mapIds(org.Hs.eg.db, keys = genes,
                  column = c('SYMBOL'), 
                  keytype = 'ENSEMBL')
genes_symbol <- data.frame(genes, symbols)

counter=1
for(gene in genes_symbol$symbols){
  if (is.na(gene)) {
    genes_symbol[counter,2] <- genes_symbol[counter,1]
  }
  counter = counter+1
}

#change ensg to gene symbols where possible
rownames(dds) <- genes_symbol$symbols
nrow(dds) 


# run DESeq2 analysis
dds <- DESeq(dds)
resultsNames(dds)  # see all comparisons 

# separate DESeq2 results by timepoints (48 hours and 10 days)
res_dif_10d <- results(dds, contrast=c("status","10d.dif","48h.undif"))
res_dif_48h <- results(dds, contrast=c("status","48h.dif","48h.undif"))

#export DESeq results in the file
#write.csv(as.data.frame(res_dif_10d[order(res_dif_10d$padj),] ), file="./ligament_dif10d_vs_undif.csv")
#write.csv(as.data.frame(res_dif_48h[order(res_dif_48h$padj),] ), file="./ligament_dif48h_vs_undif.csv")

# summary of diff.expressed genes
summary(res_dif_10d, alpha = 0.05)
summary(res_dif_48h)

# plot mean of normalized count
DESeq2::plotMA(res_dif_10d)
DESeq2::plotMA(res_dif_48h)
dev.off()

# plot counts of top gene
plotCounts(dds, gene=which.min(res_dif_10d$padj), intgroup="status")
dev.off()


## PCA
rlt <- rlog(dds)  #rlog Transformation

# tiff(file="./ligament_pca_batch.tiff",
#      units = "in",
#      width = 5,
#      height = 3,
#      res=300)
pcaData <- plotPCA(rlt, intgroup=c("status"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
names <- pcaData[,5]
pcaData[,5] <- substr(names, 1, 16)
colnames(pcaData) <-c("PC1","PC2","group","condition","replicate")
ggplot(pcaData, aes(PC1, PC2, shape=replicate, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  scale_color_brewer(palette = "Set2")
dev.off()

# remove batch effect
assay(rlt) <- limma::removeBatchEffect(assay(rlt),
                                       batch = colData(dds)[,'batch'])

# PCA after batch effect removal
# tiff(file="./ligament_pca_no_batch.tiff",
#      units = "in",
#      width = 5,
#      height = 3,
#      res=300)
pcaData <- plotPCA(rlt, intgroup=c("status"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
names <- pcaData[,5]
pcaData[,5] <- substr(names, 1, 16)
colnames(pcaData) <-c("PC1","PC2","group","condition","replicate")
ggplot(pcaData, aes(PC1, PC2, shape=replicate, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  scale_color_brewer(palette = "Set2")
dev.off()


# deg heatmap. Повышение экспрессии - красный. Понижение- синий
df <- as.data.frame(colData(dds)[,c("status")])
rownames(df) <- colnames(dds)
colnames(df) <- "status"

res_10d <- results(dds, contrast=c("status","10d.dif","48h.undif"))
# select genes with padj < 0.05
res_10d <- subset(res_10d, res_10d$padj <= 0.05 & !is.na(res_10d$padj) & abs(res_10d$log2FoldChange) >= 1.0)
# order genes by decreasing logfc
res_10d <- res_10d[order(res_10d$log2FoldChange, decreasing = TRUE),]

de <- rownames(res_10d)[1:25]  # if top-10 genes, change to rownames(res_10d)[1:10]
de_mat <- assay(rlt)[de,]
datamatrix <- t(scale(t(de_mat)))

# tiff(file="./ligament_10days_heatmap.tiff",
#      units = "in",
#      width = 6,
#      height = 7,
#      res=300)
pheatmap(datamatrix,
         show_rownames = T,
         show_colnames = T,
         annotation_col = df,
         fontsize = 7.5)
dev.off()


# volcano plot справа - усиление экспрессии, слева снижение экспрессии
res_dif10d = read.csv("./ligament_dif10d_vs_undif.csv", row.names = 1)

# tiff(file="./ligament_10dVSundif_volcano.tiff",
#      units = "in",
#      width = 8,
#      height = 6,
#      res=300)
EnhancedVolcano(res_dif10d,
                lab = rownames(res_dif10d),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                title ="MSC LIGAMENT (dif.10d vs undif)",
                labSize = 2,
                boxedLabels = F,
                col=c('black', 
                             '#CBD5E8', 
                             '#B3E2CD', 
                             '#FDCDAC'),
                             colAlpha = 1)

dev.off()


## Enrichment analysis

# combine all upregulated and downregulated genes from both timepoints
trans48 <- read.csv("./ligament_dif48h_vs_undif.csv")
trans48.up <- trans48[trans48$log2FoldChange > 1 & trans48$padj < 0.05 & !is.na(trans48$padj), 1] 
trans48.dn <- trans48[trans48$log2FoldChange < -1 & trans48$padj < 0.05 & !is.na(trans48$padj), 1]
trans48.up.entrez <- clusterProfiler::bitr(trans48.up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
trans48.dn.entrez <- clusterProfiler::bitr(trans48.dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

trans10 <- read.csv("./ligament_dif10d_vs_undif.csv")
trans10.up <- trans10[trans10$log2FoldChange > 1 & trans10$padj < 0.05 & !is.na(trans10$padj), 1] 
trans10.dn <- trans10[trans10$log2FoldChange < -1 & trans10$padj < 0.05 & !is.na(trans10$padj), 1]
trans10.up.entrez <- clusterProfiler::bitr(trans10.up,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
trans10.dn.entrez <- clusterProfiler::bitr(trans10.dn,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)


clusters <- list(trans_up_48h = trans48.up.entrez$ENTREZID,
                 trans_up_10d = trans10.up.entrez$ENTREZID,
                 trans_down_48h = trans48.dn.entrez$ENTREZID,
                 trans_down_10d = trans10.dn.entrez$ENTREZID)

#enrichGO
ck1_1 <- compareCluster(geneCluster = clusters, fun = enrichGO, OrgDb=org.Hs.eg.db, ont = "MF")

#enrichPathways

ck2_1 <- compareCluster(geneCluster = clusters, fun = enrichPathway)

#for enrichKEGG
#library(R.utils)
#R.utils::setOption("clusterProfiler.download.method","auto")

p1 <- dotplot(ck1_1, 
              title="MSC LIGAMENT upregulated, enrichGO (ONT = MF)", 
              includeAll=TRUE)+ scale_colour_distiller(palette="Set2")
p2 <- dotplot(ck2_1, 
              title="MSC LIGAMENT upregulated, enrichPathway", 
              includeAll=TRUE)+ scale_colour_distiller(palette="Set2")

# tiff(file="./ligament_ehrichGO_enrichPathway.tiff",
#      units = "in",
#      width = 15,
#      height = 20,
#      res=300)
gridExtra::grid.arrange(p1,p2,p3,p4, ncol = 2)
dev.off()


# DEG overlap between both timepoints
trans48 <- read.csv("./ligament_dif48h_vs_undif.csv")
trans48_deg <- subset(trans48, trans48$padj < 0.05 & abs(trans48$log2FoldChange) > 1)
trans48_genes <- trans48_deg$X

trans10 <- read.csv("./ligament_dif10d_vs_undif.csv")
trans10_deg <- subset(trans10, trans10$padj < 0.05 & abs(trans10$log2FoldChange) > 1)
trans10_genes <- trans10_deg$X

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(trans48_genes,
           trans10_genes),
  category.names = c("Transcriptome (48h)",
                     "Transcriptome (10d)"),
  filename = './ligament_venn.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 2000 , 
  width = 2000 , 
  resolution = 500,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  lty = 'solid',
  fill = myCol[1:2],
  
  # Numbers
  cex = .5,
  fontface = "bold",
  fontfamily = "sans",
  main.col = "black",
  
  # # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  print.mode = c("raw","percent"),
  sigdigs = 2,
  cat.fontfamily = "sans")


Reduce(intersect, list(trans48_genes, 
                       trans10_genes))
