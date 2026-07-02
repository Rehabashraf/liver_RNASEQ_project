#setup our working directory
getwd()
list.files()
#installing pacage
BiocManager::install("GEOquery")
#calliing  package  
library(tidydr)
library(dplyr)
library(GEOquery)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(IHW)
library(apeglm)
library(org.Hs.eg.db)
library(tibble)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
#load data 
livr_count=read.delim("GSE77314_raw_counts_GRCh38.p13_NCBI.tsv.gz",check.names = FALSE,row.names = 1)
livr_count=as.matrix(livr_count)
livr_meta=read.csv("SraRunTable.csv",check.names = FALSE)
#convert name of diagnosis:ch1" to condition"
colnames(livr_meta)[colnames(livr_meta)=="diagnosis"]="condition"
#match data
all(rownames(livr_meta)==colnames(livr_count))
rownames(livr_meta)=colnames(livr_count)
#convert category data to factor 
livr_meta$condition=factor(livr_meta$condition,levels = c("Normal","Tumor"))
#inspect library size
dim(livr_count)
head(livr_count[,1:4])
colSums(livr_count)/1e6
summary(colSums(livr_count)/1e6)
# visualize library size
barplot(colSums(livr_count)/1e6,las=2,ylab = "million of reads",main = "librarysize")
#log10 count distribution 
hist(log10(livr_count[livr_count>0]),breaks = 50,main = "count distribution")
#per sample
boxplot(log2(livr_count+1),las=2,outline=FALSE)
#create deseq object
type(livr_count)
dds_object=DESeqDataSetFromMatrix(countData = livr_count,colData = livr_meta,design = ~condition)
#pre filteration
table(livr_meta$condition)
smallestgroupsize=min(table(livr_meta$condition))
keep=rowSums(counts(dds_object)>=10)>=smallestgroupsize
dds_object=dds_object[keep,]
# normalization 
vsd=vst(dds_object,blind=TRUE)
#SAMPLE DISTANCE
sample_dist=dist(t(assay(vsd)))
sample_dist_matrix=as.matrix(sample_dist)
rownames(sample_dist_matrix)=paste(vsd$condition,sep="_")
colnames(sample_dist_matrix)=NULL
#visualization pheatmap
pheatmap(sample_dist_matrix,clustering_distance_rows =  sample_dist,clustering_distance_cols = sample_dist,color = colorRampPalette(rev(brewer.pal(9,"Blues")))(255),main = "samplet-to sample diatance matrix",angle_col = 45,fontsize=10)
#visualization pca
plotPCA(vsd,intgroup="condition")
pcaData<- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x=PC1, y=PC2, color=condition, label=name)) +geom_point(size=4, alpha=0.85) + geom_text(vjust=-0.8, size=3.2) +xlab(paste0("PC1:" , percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  
  theme_bw(base_size=13) + ggtitle("PCA of VST-normalized Expression") + scale_color_brewer(palette="Set1") 
# boxplot before normaliztion using condition
log_count=log2(counts(dds_object,normalized=FALSE)+1)
boxplot(log_count,las=2,col=as.integer(livr_meta$condition),main="rawlog2(counts+1)",ylab="log2(count+1)",cex.axis=.7 )
# boxplot after normaliztion using condition
dds_norm=estimateSizeFactors(dds_object)
boxplot(log2(counts(dds_norm,normalized=TRUE)+1),col=as.integer(livr_meta$condition),las=2,ylab="log2(count+1)",main="normalized count(log2+1)")
#deseq
dds=DESeq(dds_object)
sizeFactors(dds)
plotDispEsts(dds)
# extract result
livr_res=results(dds,contrast = c("condition","Tumor","Normal"),alpha = 0.05,lfcThreshold = 0,pAdjustMethod = "BH")
livr_res_df=as.data.frame(livr_res)
#extrat result using ihw
res_ihw=results(dds,contrast = c("condition","Tumor","Normal"),filterFun = ihw)
res_ihw_df=as.data.frame(res_ihw)
#shrinkage using aplegm
resultsNames(dds)
res_shrinkage=lfcShrink(dds,coef = "condition_Tumor_vs_Normal",type = "apeglm",res = livr_res)
res_shrinkage_df=as.data.frame(res_shrinkage)
# gene annotation 
res_shrinkage$genename=mapIds(org.Hs.eg.db,keys = rownames(res_shrinkage),column = "GENENAME",keytype = "ENTREZID",multiVals = "first")
res_shrinkage$ensmbl=mapIds(org.Hs.eg.db,keys = rownames(res_shrinkage),column = "ENSEMBL",keytype = "ENTREZID",multiVals = "first")
#convert to datafram for downstream analysis
res_df= as.data.frame(res_shrinkage)%>%
  rownames_to_column("gene_id")%>%
  arrange(padj)%>%
  filter(!is.na(padj))
#FILTER significant genes
sig_genes=res_df %>%
  filter(padj<0.05,abs(log2FoldChange)>1)
up_siges=filter(sig_genes,log2FoldChange>1)
down_siges=filter(sig_genes,log2FoldChange<1)
# save file
write.csv(res_df,       "DESeq2_all_results.csv",    row.names = FALSE)
write.csv(sig_genes,    "DESeq2_sig_DE_genes.csv",   row.names = FALSE)
write.csv(up_siges,     "DESeq2_upregulated.csv",    row.names = FALSE)
write.csv(down_siges,   "DESeq2_downregulated.csv",  row.names = FALSE)
#visualization MAPLOT unshinkage
plotMA(livr_res,ylim=c(-4,4), main = "MA Plot (unshrunken LFC)")
#visualization MA PLOT SHRINKAGE
plotMA(res_shrinkage,ylim=c(-4,4),main="MA Plot (shrinkage LFC)")
#volcanoplot
EnhancedVolcano(sig_genes,lab = sig_genes$ensmbl,x = "log2FoldChange",y = "padj",  title  = "Tumor vs Norml",
                subtitle       = "DESeq2 | apeglm shrinkage",
                pCutoff        = 0.05,
                FCcutoff       = 1,
                pointSize      = 2.5,
                labSize         = 3,
                col            = c("grey40", "forestgreen", "royalblue", "red2"),
                legendLabels   = c("NS", "LFC only", "p-value only", "p & LFC"),
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                max.overlaps   = 20)
# select top 50 gene for heatmape
top_genes=sig_genes %>% arrange(padj)%>%head(50)%>%
pull(gene_id)
mat=assay(vsd)[top_genes,]
mat=mat -rowMeans(mat)
ann_col <- data.frame(Condition = colData(vsd)$condition,  row.names = colnames(vsd))
pheatmap(mat,
         annotation_col   = ann_col,
         show_rownames    = TRUE,
         show_colnames    = TRUE,
         scale            = "row",
         clustering_method = "complete",
         color            = colorRampPalette(c("navy","white","firebrick3"))(100),
         fontsize_row     = 7,
         main             = "Top 50 DE Genes (DESeq2)",
         border_color     = NA
)
# enrichment 
up_entrez=up_siges$gene_id
up_entrez=na.omit(up_entrez)
# univers gene
universe=na.omit(res_df$gene_id)
# GGENNE ONTOLOOGY
ego_BP=enrichGO(gene = up_entrez,universe = universe,OrgDb =org.Hs.eg.db,ont = "BP",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,readable = TRUE )
ego_MF=enrichGO(gene = up_entrez,universe = universe,OrgDb =org.Hs.eg.db,ont = "MF",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,readable = TRUE )
ego_CC=enrichGO(gene = up_entrez,universe = universe,OrgDb =org.Hs.eg.db,ont = "CC",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,readable = TRUE )
ego_ALL=enrichGO(gene = up_entrez,universe = universe,OrgDb =org.Hs.eg.db,ont = "ALL",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.2,readable = TRUE )
# SAVEE FILEE
write.csv(ego_BP,"biological.csv",row.names = FALSE)
write.csv(ego_MF,"molecular.csv",row.names = FALSE)
write.csv(ego_CC,"cellular.csv",row.names = FALSE)
write.csv(ego_BP,"all.csv",row.names = FALSE)
# visualiztion dotplot
dotplot(ego_BP, showCategory = 20, title = "GO Biological Process")
dotplot(ego_MF, showCategory = 20, title = "GO Molecular function")
dotplot(ego_CC, showCategory = 15, title = "GO cellular component")
# visualiztion barplot
barplot(ego_BP, showCategory = 10)
barplot(ego_MF, showCategory = 10)
barplot(ego_CC, showCategory = 10)
# visualiztion emaplot up
emapplot(pairwise_termsim(ego_BP),showCategory=10)
emapplot(pairwise_termsim(ego_MF),showCategory=10)
emapplot(pairwise_termsim(ego_CC),showCategory=10)
# kgg pathway
kggpathway_up=enrichKEGG(gene = up_entrez,organism = "hsa",pvalueCutoff = 0.05,universe = universe)
#save kggpathway
write.csv(kggpathway_up,"kggpathway_up.csv",row.names = FALSE)
#visualization kggpathway
dotplot(kggpathway_up,showCategory=10,title="kggpathway_up")
barplot(kggpathway_up, showCategory = 10)
emapplot(pairwise_termsim(kggpathway_up),showCategory=10)
##GSEA
ranked_gene=sig_genes%>%
  filter(!is.na(gene_id),!is.na(pvalue))%>%
  mutate(rank=-log10(pvalue)*sign(log2FoldChange))%>%
  arrange(desc(rank))
ranked_gene_list=ranked_gene$rank
names(ranked_gene_list)=ranked_gene$gene_id
gsea_res <- gseKEGG(
  geneList      = ranked_gene_list,
  organism      = "hsa",
  nPermSimple   = 1000,
  minGSSize     = 15,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  verbose       = FALSE)

#visualization GSEA
gseaplot(gsea_res,geneSetID = 1,title = gsea_res$Description[1])
dotplot(gsea_res)
ridgeplot(gsea_res,showCategory = 12)
emapplot(pairwise_termsim(gsea_res),showCategory=10)
#Reactome up__gene
reactom_res_up=enrichPathway(
  gene     = up_entrez,
  organism = "human",
  universe =universe,
  pvalueCutoff = 0.05)
#save file 
write.csv(reactom_res_up,"reactom_res_up.csv",row.names = FALSE)
#visuualizatiion for  up  
dotplot(reactom_res_up,showCategory=10,title="reactom_res_up")
barplot(reactom_res_up, showCategory = 10)
emapplot(pairwise_termsim(reactom_res_up),showCategory=10)


