library(rtracklayer)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(reshape2)

p53 <- read.delim("p53_signaling_pathway.txt", skip = 2, header = F)
mismatch <- read.delim("Mismatch_repair.txt", skip = 2, header = F)
brca.counts <- read.csv("countdata.csv")
brca.counts2 <- brca.counts[,-1]
rownames(brca.counts2) <- brca.counts[,1]

## Annotate gene names
gtf = readGFF("gencode.v37.annotation.gtf")
gtf = gtf[gtf$type == 'gene',]
gtf[gtf$gene_name %in% p53$V1,]$gene_id

p53.all <- c(p53$V1, mismatch$V1)
brca.path <- brca.counts2[rownames(brca.counts2) %in% gtf[gtf$gene_name %in% p53.all,]$gene_name, ]

p53.brca <- brca.counts2[rownames(brca.counts2) %in% gtf[gtf$gene_name %in% p53$V1, ]$gene_name, ]

mismatch.brca <- brca.counts2[rownames(brca.counts2) %in% gtf[gtf$gene_name %in% mismatch$V1, ]$gene_name, ]


## Heatmap
metadata <- data.frame(
  factor(c(rep("FT282_Control",4), rep("FT282_BRCA1_OCR",4), rep("FT282_BRCA2_BCR",2),
           rep("FT282_BRCA2_OCR",4), rep("MCF10A_BRCA1_BCR",2), rep("MCF10A_BRCA2_BCR",4),
           rep("MCF10A_Control",4), rep("MCF10A_BRCA1_OCR",4))),
  row.names = colnames(p53.brca))
colnames(metadata) <- "condition"
  
pheatmap(p53.brca %>% t %>% scale %>% t %>% na.omit(), , border_color = NA, color = viridis(50), show_rownames = F,
         main = "p53 signaling pathway genes", annotation_col = metadata)

pheatmap(mismatch.brca %>% t %>% scale %>% t %>% na.omit(), , border_color = NA, color = viridis(50), show_rownames = F,
         main = "mismatch repair pathway genes", annotation_col = metadata)

## Boxplot of different pathway genes
box <- as.data.frame(t(p53.brca))
Group <- c(rep("FT282",14), rep("MCF10A",14))
brca.anno <- data.frame(Group)
box$Group <- brca.anno$Group


m.box <- melt(box)
ggplot(m.box, aes(Group, log1p(value), fill = Group)) + geom_boxplot() + 
  scale_fill_manual(values = brewer.pal(2, "Set2")) + 
  ggtitle("P53 signaling genes expression") + stat_compare_means()

box2 <- as.data.frame(t(mismatch.brca))
box2$Group <- brca.anno$Group

m.box2 <- melt(box2)
ggplot(m.box2, aes(Group, log1p(value), fill = Group)) + geom_boxplot() + 
  scale_fill_manual(values = brewer.pal(2, "Set2")) + 
  ggtitle("mismatch repair genes expression") + stat_compare_means()

## Boxplot for conditions in p53 signalling pathway
box.1 <- box[c(5:8,1:4),]
Group.1 <- c(rep("FT282-BRCA1_OCR",4),rep("Control",4)) 
box.1_anno <- data.frame(Group.1)
box.1$Group <- box.1_anno$Group.1

box.2 <- box[c(11:14,1:4),]
Group.2 <- c(rep("FT282-BRCA2_OCR",4),rep("Control",4))
box.2_anno <- data.frame(Group.2)
box.2$Group <- box.2_anno$Group.2

box.3 <- box[c(5:8,11:14),]
Group.3 <- c(rep("FT282-BRCA1_OCR",4),rep("FT282-BRCA2_OCR",4))
box.3_anno <- data.frame(Group.3)
box.3$Group <- box.3_anno$Group.3

pbox.1 <- melt(box.1)
pbox.2 <- melt(box.2)
pbox.3 <- melt(box.3)

bp1 <- ggplot(pbox.1, aes(Group, log1p(value), fill = Group)) + geom_boxplot() + 
  scale_fill_manual(values = brewer.pal(2, "Set2")) + 
  ggtitle(expression("FT282"^{BRCA1:OCR}*"mutations vs Controls")) + stat_compare_means()

bp2 <- ggplot(pbox.2, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("FT282"^{BRCA2:OCR}*"mutations vs Controls")) + stat_compare_means()

bp3 <- ggplot(pbox.3, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("FT282"^{BRCA1:OCR}*"mutations vs FT282"^{BRCA2:OCR}*"mutations")) + stat_compare_means()

## Multiple plots same page
fig1 <- ggarrange(bp1, bp2, bp3, labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
annotate_figure(fig1, top = text_grob("P53 signaling pathway", color = "red" ,face = "bold"))

box.4 <- box[c(25:28,21:24),]
Group.4 <- c(rep("MCF10A-BRCA1_OCR",4),rep("Control",4)) 
box.4_anno <- data.frame(Group.4)
box.4$Group <- box.4_anno$Group.4

box.5 <- box[c(17:20,21:24),]
Group.5 <- c(rep("MCF10A-BRCA2_BCR",4),rep("Control",4)) 
box.5_anno <- data.frame(Group.5)
box.5$Group <- box.5_anno$Group.5

box.6 <- box[c(25:28,17:20),]
Group.6 <- c(rep("MCF10A-BRCA1_OCR",4),rep("MCF10A-BRCA2_BCR",4)) 
box.6_anno <- data.frame(Group.6)
box.6$Group <- box.6_anno$Group.6

pbox.4 <- melt(box.4)
pbox.5 <- melt(box.5)
pbox.6 <- melt(box.6)

bp4 <- ggplot(pbox.4, aes(Group, log1p(value), fill = Group)) + geom_boxplot() + 
  scale_fill_manual(values = brewer.pal(2, "Set2")) + 
  ggtitle(expression("MCF10A"^{BRCA1:OCR}*"mutations vs Controls")) + stat_compare_means()

bp5 <- ggplot(pbox.5, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("MCF10A"^{BRCA2:BCR}*"mutations vs Controls")) + stat_compare_means()

bp6 <- ggplot(pbox.6, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("MCF10A"^{BRCA1:OCR}*"mutations vs MCF10A"^{BRCA2:BCR}*"mutations")) + stat_compare_means()

## Multiple plots same page
fig2 <- ggarrange(bp4, bp5, bp6, labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
annotate_figure(fig2, top = text_grob("P53 signaling pathway", color = "red" ,face = "bold"))


box.7 <- box[c(5:8,25:28),]
Group.7 <- c(rep("FT282-BRCA1_OCR",4),rep("MCF10A-BRCA1_OCR",4)) 
box.7_anno <- data.frame(Group.7)
box.7$Group <- box.7_anno$Group.7

box.8 <- box[c(11:14,17:20),]
Group.8 <- c(rep("FT282-BRCA2_OCR",4),rep("MCF10A-BRCA2_BCR",4)) 
box.8_anno <- data.frame(Group.8)
box.8$Group <- box.8_anno$Group.8

pbox.7 <- melt(box.7)
pbox.8 <- melt(box.8)

bp7 <- ggplot(pbox.7, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("FT282"^{BRCA1:OCR}*"mutations vs MCF10A"^{BRCA1:OCR}*"mutations")) + stat_compare_means()

bp8 <- ggplot(pbox.8, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("FT282"^{BRCA2:OCR}*"mutations vs MCF10A"^{BRCA2:BCR}*"mutations")) + stat_compare_means()

## Multiple plots same page
fig3 <- ggarrange(bp7, bp8, labels = c("A", "B"),
          ncol = 2, nrow = 2)
annotate_figure(fig3, top = text_grob("P53 signaling pathway", color = "red" ,face = "bold"))

## Boxplot for conditions in Mismatch repair pathway
box2$Group <- NULL
box.9 <- box2[c(5:8,1:4),]
Group.9 <- c(rep("FT282-BRCA1_OCR",4),rep("Control",4)) 
box.9_anno <- data.frame(Group.9)
box.9$Group <- box.9_anno$Group.9

box.10 <- box2[c(11:14,1:4),]
Group.10 <- c(rep("FT282-BRCA2_OCR",4),rep("Control",4)) 
box.10_anno <- data.frame(Group.10)
box.10$Group <- box.10_anno$Group.10

box.11 <- box2[c(5:8,11:14),]
Group.11 <- c(rep("FT282-BRCA1_OCR",4),rep("FT282-BRCA2_OCR",4)) 
box.11_anno <- data.frame(Group.11)
box.11$Group <- box.11_anno$Group.11

pbox.9 <- melt(box.9)
pbox.10 <- melt(box.10)
pbox.11 <- melt(box.11)

bp9 <- ggplot(pbox.9, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("FT282"^{BRCA1:OCR}*"mutations vs Controls")) + stat_compare_means()

bp10 <- ggplot(pbox.10, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("FT282"^{BRCA2:OCR}*"mutations vs Controls")) + stat_compare_means()

bp11 <- ggplot(pbox.11, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("FT282"^{BRCA1:OCR}*"mutations vs FT282"^{BRCA2:OCR}*"mutations")) + stat_compare_means()

## Multiple plots same page
fig4 <- ggarrange(bp9, bp10, bp11, labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
annotate_figure(fig4, top = text_grob("Mismatch repair pathway", color = "red" ,face = "bold"))


box.12 <- box2[c(25:28,21:24),]
Group.12 <- c(rep("MCF10A-BRCA1_OCR",4),rep("Control",4)) 
box.12_anno <- data.frame(Group.12)
box.12$Group <- box.12_anno$Group.12

box.13 <- box2[c(17:20,21:24),]
Group.13 <- c(rep("MCF10A-BRCA2_BCR",4),rep("Control",4)) 
box.13_anno <- data.frame(Group.13)
box.13$Group <- box.13_anno$Group.13

box.14 <- box2[c(25:28,17:20),]
Group.14 <- c(rep("MCF10A-BRCA1_OCR",4),rep("MCF10A-BRCA2_BCR",4)) 
box.14_anno <- data.frame(Group.14)
box.14$Group <- box.14_anno$Group.14

pbox.12 <- melt(box.12)
pbox.13 <- melt(box.13)
pbox.14 <- melt(box.14)

bp12 <- ggplot(pbox.12, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("MCF10A"^{BRCA1:OCR}*"mutations vs Controls")) + stat_compare_means()

bp13 <- ggplot(pbox.13, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("MCF10A"^{BRCA2:BCR}*"mutations vs Controls")) + stat_compare_means()

bp14 <- ggplot(pbox.14, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("MCF10A"^{BRCA1:OCR}*"mutations vs MCF10A"^{BRCA2:BCR}*"mutations")) + stat_compare_means()

## Multiple plots same page
fig5 <- ggarrange(bp12, bp13, bp14, labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)
annotate_figure(fig5, top = text_grob("Mismatch repair pathway", color = "red" ,face = "bold"))


box.15 <- box2[c(5:8,25:28),]
Group.15 <- c(rep("FT282-BRCA1_OCR",4),rep("MCF10A-BRCA1_OCR",4)) 
box.15_anno <- data.frame(Group.15)
box.15$Group <- box.15_anno$Group.15

box.16 <- box2[c(11:14,17:20),]
Group.16 <- c(rep("FT282-BRCA2_OCR",4),rep("MCF10A-BRCA2_BCR",4)) 
box.16_anno <- data.frame(Group.16)
box.16$Group <- box.16_anno$Group.16

pbox.15 <- melt(box.15)
pbox.16 <- melt(box.16)

bp15 <- ggplot(pbox.15, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("FT282"^{BRCA1:OCR}*"mutations vs MCF10A"^{BRCA1:OCR}*"mutations")) + stat_compare_means()

bp16 <- ggplot(pbox.16, aes(Group, log1p(value), fill = Group)) + geom_boxplot() +
  scale_fill_manual(values = brewer.pal(2, "Set2")) +
  ggtitle(expression("FT282"^{BRCA2:OCR}*"mutations vs MCF10A"^{BRCA2:BCR}*"mutations")) + stat_compare_means()

## Multiple plots same page
fig6 <- ggarrange(bp15, bp16, labels = c("A", "B"),
          ncol = 2, nrow = 2)
annotate_figure(fig6, top = text_grob("Mismatch repair pathway", color = "red" ,face = "bold"))
