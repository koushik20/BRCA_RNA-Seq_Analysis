# Loading Packages
library(dplyr)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(ggrepel)
library(viridis)

countdata <- read.csv("countdata.csv")
countdata <- data.frame(countdata[,-1], row.names=countdata[,1])
countdata <- as.matrix(countdata)
head(countdata)

# new dataframe for downstream analysis
countdata2<- countdata[ , c(1:4,11:14)]
dim(countdata2)

# Assign condition (first four are controls, second four contain the expansion)
(condition <- factor(c(rep("FT282_OR10A4",4), rep("FT282_BRCA2OCR", 4))))

# Analysis with DESeq2
(coldata <- data.frame(row.names=colnames(countdata2), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata2, colData=coldata, design=~condition)
dds

# Running DESeq pipeline
dds <- DESeq(dds)
res <- results(dds)
summary(res)
res

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
hist(assay(rld))

# Principal components analysis
png("qc-pca.png", 1000, 1000, pointsize=20)
DESeq2::plotPCA(rld, intgroup="condition")
dev.off()

################################################################################
################################################################################

# Heatmap of the 100 most variant genes
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 100)
matrix <- assay(rld)[ topVarGenes, ]
matrix <- matrix - rowMeans(matrix)
anno <- as.data.frame(colData(rld))
pheatmap(matrix, annotation_col = anno, scale = "row", border_color = NA, fontsize = 6, main = "Control vs BRCA2 OCR1")

################################################################################
################################################################################

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

################################################################################
################################################################################


# HeatMap of DEG's
de <- rownames(res[res$padj<0.05 & !is.na(res$padj), ])
de_mat <- assay(rld)[de, ] # Differentially expressed gene list
anno <- as.data.frame(colData(rld))
anno2 <- subset(anno, select = condition)
pheatmap(de_mat,show_rownames = F,  color = viridis(15, direction = -1), 
         annotation_col = anno2, scale = "row", border_color = NA, fontsize = 6, 
         main = expression("FT282"^{BRCA2:OCR}*"mutations vs Controls"))


################################################################################
################################################################################

# Volcano Plot
res2 <- na.omit(res)
res2$STATUS <- "Not Sig"
res2[res2$log2FoldChange < -1 & res2$padj < 0.05, ]$STATUS <- "Down"
res2[res2$log2FoldChange > 1 & res2$padj < 0.05, ]$STATUS <- "Up"

aux <- as.data.frame(res2)
aux <- cbind(Gene_ID = rownames(aux), aux)
rownames(aux) <- 1:nrow(aux)
aux$Symbol <- aux$Gene_ID

# label top DEGs on volcano
aux$label <- NA
up <- aux[aux$STATUS == "Up", ]
up <- arrange(up, -log2FoldChange, padj)
up <- arrange(up,  padj)

up$label[1:10] <- up$Symbol[1:10]

dn <- aux[aux$STATUS == "Down", ]
dn <- arrange(dn, log2FoldChange, padj)
dn <- arrange(dn, padj)

dn$label[1:10] <- dn$Symbol[1:10]

# transfer label to the plotting obj
aux$label[which(aux$Symbol %in% up$label[1:10])] <- up$label[1:10]
aux$label[which(aux$Symbol %in% dn$label[1:10])] <- dn$label[1:10]


pdf("Volcano_BRCA2_OCRvsControl.pdf", 15, 10)
ggplot(aux, aes(log2FoldChange, -log10(padj), col = STATUS, label = label)) + geom_point() + 
  geom_label_repel(aes(label = aux$label), vjust = "inward", hjust = "inward") + xlab("log2 Fold Change") +
  ylab("-log10 Adjusted Pvalue") + ggtitle(expression("FT282"^{BRCA2:OCR}*"mutations vs Controls")) +theme_bw() +
  scale_color_manual(values = c("darkolivegreen4",  "grey", "deepskyblue4"))
dev.off()


################################################################################
################################################################################

# MA Plot
# Get differential expression results
res <- results(dds)
table(res$padj<0.05)

# Order by adjusted p-value
res <- res[order(res$padj), ]

# Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

# Write results
write.csv(aux, file="diffexpr-results_Control_vs_BRCA2.csv")

# Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

# Examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

# MA plot
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot-2.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot between Control vs BRCA2 OCR1")
dev.off()

