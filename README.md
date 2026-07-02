# Liver Cancer RNA-Seq Analysis Using DESeq2

## Objectives 


Liver cancer is one of the most common causes of cancer-related deaths worldwide. Understanding gene expression changes between tumor and normal tissues can provide valuable insights into the molecular mechanisms of hepatocellular carcinoma (HCC).

This project aims to :

- Identify differentially expressed genes (DEGs) between tumor and normal liver tissues.

- Investigate biological processes associated with liver cancer progression.

- Explore enriched molecular pathways and signaling networks.



## Dataset

- Dataset ID:GSE77314

-Source: GEO 

### Study Design

Samples include:

- Normal liver tissue

- Liver tumor tissue

For differential expression analysis, the comparison performed in this project is:

- Tumor vs Normal

## Bioinformatics Workflow 

### Data Import


- RNA-Seq count matrix loading

- Sample metadata loading

- Metadata matching with count matrix


### Data Preprocessing

- Conversion of condition to factors

- Construction of DESeq2 dataset

### Gene Filtering


Lowly expressed genes were removed using:

- Minimum count threshold ≥ 10

- Expression required in at least the smallest experimental group



### Quality Control

- Library size inspection

- Count distribution visualization

- Boxplots before normalization

- Sample-to-sample distance heatmap

- Principal Component Analysis (PCA)


### Normalization

Variance Stabilizing Transformation VST


### Differential Expression Analysis

Differential expression analysis was performed using DESeq2.

Comparison:

Tumor vs Normal

Statistical criteria:

- Adjusted p-value < 0.05

- |log2 Fold Change| > 1


### Gene Annotation

Gene identifiers were mapped to:

- Gene Names

- Ensembl IDs

using org.Hs.eg.db


### Visualization

Generated visualizations include:

- PCA plots

- Sample distance heatmaps

- MA plots

- Volcano plots

- Heatmaps of differentially expressed genes

### Functional Enrichment Analysis

Gene Ontology (GO):

- Biological Process (BP)

- Molecular Function (MF)

- Cellular Component (CC)


### Pathway Analysis

KEGG pathway enrichment analysis was performed to identify significantly enriched biological pathways.


### Gene Set Enrichment Analysis (GSEA)

Ranked-gene GSEA was performed to identify coordinated pathway-level changes associated with liver cancer.


### Reactome Analysis

Reactome pathway enrichment analysis was performed to investigate higher-order biological mechanisms.


###  Packages

- DESeq2

- ggplot2

- pheatmap

- EnhancedVolcano

- clusterProfiler

- ReactomePA

- enrichplot

- org.Hs.eg.db

- dplyr

- apeglm

- IHW

### Programming Language

- R



## Output Files

### Differential Expression Results


- DESeq2_all_results.csv

- DESeq2_sig_DE_genes.csv

- DESeq2_upregulated.csv

- DESeq2_downregulated.csv


### GO Enrichment

- biological.csv

- molecular.csv

- cellular.csv


### KEGG Analysis


- kggpathway_up.csv


### Reactome Analysis


- reactom_res_up.csv















