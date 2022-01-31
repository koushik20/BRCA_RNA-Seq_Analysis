
## Pathway Analysis of all DEGs Control vs BRCA1 OCR1 
library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

res3 <- read.csv("diffexpr-results_Control_vs_BRCA2.csv")
res3 <- data.frame(res3[,-1], row.names=res3[,1])
res3
norm.can.path <- res3[res3$padj < 0.05, ]$log2FoldChange

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes.mart <- 
  getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
        values = substr(rownames(res3[res3$padj < 0.05, ]),1,15), mart = mart)

res3$entrezid <- genes.mart$entrezgene_id[match(substr(rownames(res3),1,15), genes.mart$ensembl_gene_id)]
names(norm.can.path) <- res3[res3$padj < 0.05,]$entrezid


gene.df <- bitr(names(norm.can.path), fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
gene.df <- gene.df[gene.df$ENSEMBL %in% substr(rownames(res3),1,15),]
aux2 <- res3[substr(rownames(res3),1,15) %in% gene.df$ENSEMBL, ]
gene.df$FC <- aux2$log2FoldChange[match(substr(rownames(aux2),1,15), gene.df$ENSEMBL)]

ego <- enrichGO(gene         = gene.df$ENSEMBL,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
head(ego)


ego.mf <-  enrichGO(gene         = gene.df$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)


ego.cc <-  enrichGO(gene         = gene.df$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)


e.go.bp <- setReadable(ego, 'org.Hs.eg.db', 'ENSEMBL')
e.go.mf <- setReadable(ego.mf, 'org.Hs.eg.db', 'ENSEMBL')
e.go.cc <- setReadable(ego.cc, 'org.Hs.eg.db', 'ENSEMBL')


write.table(e.go.bp, quote =F, row.names = F, sep = "\t", 
            file = "ControlvsBRCA2_OCR1_GO_BP_table.tsv")
write.table(e.go.mf, quote =F, row.names = F, sep = "\t", 
            file = "ControlvsBRCA2_OCR1_GO_MF_table.tsv")
write.table(e.go.cc, quote =F, row.names = F, sep = "\t", 
            file = "ControlvsBRCA2_OCR1_GO_CC_table.tsv")

#DotPlots
library(enrichplot)
pdf("Pathway_DotPlot_ControlvsBRCA2_OCR_GOBP.pdf", 17, 15)
dotplot(ego, showCategory=30, title = expression("FT282"^{BRCA2:OCR}*"mutations vs Controls"))
dev.off()

pdf("Pathway_Dotplot_ControlvsBRCA2_OCR_GOMF.pdf", 17, 15)
dotplot(ego.mf, showCategory=50,  title = expression("FT282"^{BRCA2:OCR}*"mutations vs Controls"))
dev.off()

pdf("Pathway_DotPlot_ControlvsBRCA2_OCR1_GOCC.pdf", 17, 15)
dotplot(ego.cc, showCategory = 50, title = "Control vs BRCA2 OCR1 , GO CC")
dev.off()

#Goplots
pdf("Pathway_GOPlot_ControlvsBRCA2_OCR1_GOB.pdf", 20, 15)
goplot(ego, showCategory=20, title = "Control vs BRCA2 OCR1 , GO BP")
dev.off()

pdf("Pathway_GOPlot_ControlvsBRCA2_OCR1_GOMF.pdf", 20, 15)
goplot(ego.mf, showCategory=20, title = "Control vs BRCA2 OCR1 , GO MF")
dev.off()

pdf("Pathway_GOPlot_ControlvsBRCA2_OCR1_GOCC.pdf", 20, 15)
goplot(ego.cc, showCategory=20, title = "Control vs BRCA2 OCR1 , GO CC")
dev.off()

################################################################################
################################################################################

# Pathway analysis for upregulated genes of Control vs BRCA2 OCR1
library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

res3 <- read.csv("res_data_ControlvsBRCA2.csv")
res3 <- data.frame(res3[,-1], row.names=res3[,1])
res3

library(dplyr)
res4 <- filter(res3, STATUS == "Up")

norm.can.path <- res4[res4$padj < 0.05, ]$log2FoldChange

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes.mart <- 
  getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
        values = substr(rownames(res4[res4$padj < 0.05, ]),1,15), mart = mart)

res4$entrezid <- genes.mart$entrezgene_id[match(substr(rownames(res4),1,15), genes.mart$ensembl_gene_id)]
names(norm.can.path) <- res4[res4$padj < 0.05,]$entrezid


gene.df <- bitr(names(norm.can.path), fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
gene.df <- gene.df[gene.df$ENSEMBL %in% substr(rownames(res4),1,15),]
aux3 <- res4[substr(rownames(res4),1,15) %in% gene.df$ENSEMBL, ]
gene.df$FC <- aux3$log2FoldChange[match(substr(rownames(aux3),1,15), gene.df$ENSEMBL)]

ego <- enrichGO(gene         = gene.df$ENSEMBL,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
head(ego)


ego.mf <-  enrichGO(gene         = gene.df$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)


ego.cc <-  enrichGO(gene         = gene.df$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)


e.go.bp <- setReadable(ego, 'org.Hs.eg.db', 'ENSEMBL')
e.go.mf <- setReadable(ego.mf, 'org.Hs.eg.db', 'ENSEMBL')
e.go.cc <- setReadable(ego.cc, 'org.Hs.eg.db', 'ENSEMBL')


write.table(e.go.bp, quote =F, row.names = F, sep = "\t", 
            file = "ControlvsBRCA2_OCR1_GO_BP_table_Upregulated.tsv")
write.table(e.go.mf, quote =F, row.names = F, sep = "\t", 
            file = "ControlvsBRCA2_OCR1_GO_MF_table_Upregulated.tsv")
write.table(e.go.cc, quote =F, row.names = F, sep = "\t", 
            file = "ControlvsBRCA2_OCR1_GO_CC_table_Upregulated.tsv")

#DotPlots
library(enrichplot)
pdf("Pathway_DotPlot_ControlvsBRCA2_OCR1_GOBP_Upregulated.pdf", 17, 15)
dotplot(ego, showCategory=50, title = "Control vs BRCA2 OCR1 of Upregulated Genes, GO BP")
dev.off()

pdf("Pathway_Dotplot_ControlvsBRCA2_OCR1_GOMF_Upregulated.pdf", 17, 15)
dotplot(ego.mf, showCategory=50,  title = "Control vs BRCA2 OCR1 of Upregulated Genes, GO MF")
dev.off()

pdf("Pathway_DotPlot_ControlvsBRCA2_OCR1_GOCC_Upregulated.pdf", 17, 15)
dotplot(ego.cc, showCategory = 50, title = "Control vs BRCA2 OCR1 of Upregulated Genes, GO CC")
dev.off()

#Goplots
pdf("Pathway_GOPlot_ControlvsBRCA2_OCR1_GOBP_Upregulated.pdf", 20, 15)
goplot(ego, showCategory=20, title = "Control vs BRCA2 OCR1 of Upregulated Genes, GO BP")
dev.off()

pdf("Pathway_GOPlot_ControlvsBRCA2_OCR1_GOMF_Upregulated.pdf", 20, 15)
goplot(ego.mf, showCategory=20, title = "Control vs BRCA2 OCR1 of Upregulated Genes, GO MF")
dev.off()

pdf("Pathway_GOPlot_ControlvsBRCA2_OCR1_GOCC_Upregulated.pdf", 20, 15)
goplot(ego.cc, showCategory=20, title = "Control vs BRCA2 OCR1 of Upregulated Genes, GO CC")
dev.off()


################################################################################
################################################################################

# Pathway analysis for downregulated genes of Control vs BRCA2 OCR1
library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

res3 <- read.csv("res_data_ControlvsBRCA2.csv")
res3 <- data.frame(res3[,-1], row.names=res3[,1])
res3

library(dplyr)
res5 <- filter(res3, STATUS == "Down")

norm.can.path <- res5[res5$padj < 0.05, ]$log2FoldChange

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes.mart <- 
  getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
        values = substr(rownames(res5[res5$padj < 0.05, ]),1,15), mart = mart)

res5$entrezid <- genes.mart$entrezgene_id[match(substr(rownames(res5),1,15), genes.mart$ensembl_gene_id)]
names(norm.can.path) <- res5[res5$padj < 0.05,]$entrezid


gene.df <- bitr(names(norm.can.path), fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
gene.df <- gene.df[gene.df$ENSEMBL %in% substr(rownames(res5),1,15),]
aux4 <- res5[substr(rownames(res5),1,15) %in% gene.df$ENSEMBL, ]
gene.df$FC <- aux4$log2FoldChange[match(substr(rownames(aux4),1,15), gene.df$ENSEMBL)]

ego <- enrichGO(gene         = gene.df$ENSEMBL,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
head(ego)


ego.mf <-  enrichGO(gene         = gene.df$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)


ego.cc <-  enrichGO(gene         = gene.df$ENSEMBL,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)


e.go.bp <- setReadable(ego, 'org.Hs.eg.db', 'ENSEMBL')
e.go.mf <- setReadable(ego.mf, 'org.Hs.eg.db', 'ENSEMBL')
e.go.cc <- setReadable(ego.cc, 'org.Hs.eg.db', 'ENSEMBL')


write.table(e.go.bp, quote =F, row.names = F, sep = "\t", 
            file = "ControlvsBRCA2_OCR1_GO_BP_table_downregulated.tsv")
write.table(e.go.mf, quote =F, row.names = F, sep = "\t", 
            file = "ControlvsBRCA2_OCR1_GO_MF_table_downregulated.tsv")
write.table(e.go.cc, quote =F, row.names = F, sep = "\t", 
            file = "ControlvsBRCA2_OCR1_GO_CC_table_downregulated.tsv")

#DotPlots
library(enrichplot)
pdf("Pathway_DotPlot_ControlvsBRCA2_OCR1_GOBP_downregulated.pdf", 17, 15)
dotplot(ego, showCategory=50, title = "Control vs BRCA2 OCR1 of downregulated Genes, GO BP")
dev.off()

pdf("Pathway_Dotplot_ControlvsBRCA2_OCR1_GOMF_downregulated.pdf", 17, 15)
dotplot(ego.mf, showCategory=50,  title = "Control vs BRCA2 OCR1 of downregulated Genes, GO MF")
dev.off()

pdf("Pathway_DotPlot_ControlvsBRCA2_OCR1_GOCC_downregulated.pdf", 17, 15)
dotplot(ego.cc, showCategory = 50, title = "Control vs BRCA2 OCR1 of downregulated Genes, GO CC")
dev.off()

#Goplots
pdf("Pathway_GOPlot_ControlvsBRCA2_OCR1_GOBP_downregulated.pdf", 20, 15)
goplot(ego, showCategory=20, title = "Control vs BRCA2 OCR1 of downregulated Genes, GO BP")
dev.off()

pdf("Pathway_GOPlot_ControlvsBRCA2_OCR1_GOMF_downregulated.pdf", 20, 15)
goplot(ego.mf, showCategory=20, title = "Control vs BRCA2 OCR1 of downregulated Genes, GO MF")
dev.off()

pdf("Pathway_GOPlot_ControlvsBRCA2_OCR1_GOCC_downregulated.pdf", 20, 15)
goplot(ego.cc, showCategory=20, title = "Control vs BRCA2 OCR1 of downregulated Genes, GO CC")
dev.off()

