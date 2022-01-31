# BRCA_RNA-Seq_Analysis
Code repository for the manuscript "Cell type specific functional effects of different BRCA1 and BRCA2 pathogenic mutations in breast and ovarian cancer precursor cells"

Differential Expression Analysis of 8 different conditions

1. FT282-BRCA1_OCR vs Control
2. FT282-BRCA2_OCR vs Control
3. FT282-BRCA1_OCR vs FT282-BRCA2_OCR
4. MCF10A-BRCA1_OCR vs Control
5. MCF10A-BRCA2_BCR vs Control
6. MCF10A-BRCA1_OCR vs MCF10A-BRCA2_BCR
7. MCF10A-BRCA1_OCR vs FT282-BRCA1_OCR
8. MCF10A-BRCA2_BCR vs FT282-BRCA2_OCR

#### Files ####

**Countdata.csv** - Raw counts used in the statistical programs downstream for differential gene expression where each row is a gene and eash column is a sample.

**DESeq2_FT282-BRCA2_OCR vs Controls.r** - R code for Differential expression analysis using countdata.

**diffexpr-results_Control_vs_BRCA2.csv** - Differentially expressed genes results generated from DESeq2 and further used for enrichment analysis.

**Enrichment_Analysis.R** - Enrihcemnt analysis that generates DotPlots and GOPlots with enrichment GO categories after FDR control.

**Mismatch_repair.txt** - Mismatch repair pathway genes from KEGG database.

**p53_signaling_pathway.txt** - p53 signaling pathway genes from KEGG database.

**pathway_analysis.R** - Specific pathway analysis 

* We performed differential expression analysis for each of the eight conditions using the same workflow.

