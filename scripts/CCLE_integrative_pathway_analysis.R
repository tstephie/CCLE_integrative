# CCLE Integrative Pathway Analysis
## RNA-seq and RPPA
## 3/11/20

library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)
library(ggplot2)
library(DESeq2)

# set working directory
setwd('~/tutorial_r/integrative/CCLE')

# read in data 
ccle_integrate <- read.csv('final_CCLE_analysis_v2.csv', header = T, stringsAsFactors = F)
ccle_integrate <- ccle_integrate[,-1]

# check summaries
summary(ccle_integrate$max_log_FC)

# Get pathways database
m_df <- msigdbr(species = 'Homo sapiens', category = 'H') %>% dplyr::select(gs_name, entrez_gene)
m_df_filt <- m_df %>% dplyr::select(gs_name, entrez_gene)
m_t2g <- msigdbr('Homo sapiens', category = 'C2') %>% dplyr::select(gs_name, entrez_gene) # C2 pathways

# Convert genes symbols into Entrez IDs
gene_convert <- bitr(ccle_integrate$gene_name, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Hs.eg.db", drop = F)

# filter dataset with matching Entrez IDS
ccle_integrate <- ccle_integrate %>% mutate(ENTREZID = gene_convert$ENTREZID) %>% filter(!is.na(ENTREZID))

# make named vector of FC 
gene_list <- ccle_integrate$max_log_FC
names(gene_list) <- ccle_integrate$ENTREZID
gene_list <- sort(gene_list, decreasing = T)

# Run GSEA
ccle_integrate_gsea_C2 <- GSEA(gene_list, TERM2GENE = m_df, pvalueCutoff = 1)
ccle_integrate_gsea_C2@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), NES, fill = p.adjust)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("Integrative GSEA Pathways") + xlab(element_blank()) + ylab('NES') + scale_fill_gradient(low = 'blue', high = 'red')

# # Run GSEA with filtered FC
# ccle_integrate_filt <- ccle_integrate %>% dplyr::filter(abs(max_log_FC) > 2)
# gene_list <- ccle_integrate_filt$max_log_FC
# names(gene_list) <- ccle_integrate_filt$ENTREZID
# gene_list <- sort(gene_list, decreasing = T)
# ccle_integrate_filt_gsea_C2 <- GSEA(gene_list, TERM2GENE = m_t2g, pvalueCutoff = 1)

# Run GSEA KEGG
gene_list <- ccle_integrate$max_log_FC
names(gene_list) <- ccle_integrate$ENTREZID
gene_list <- sort(gene_list, decreasing = T)
ccle_integrate_kegg <- gseKEGG(gene_list, pvalueCutoff = 1)
## description plot
ccle_integrate_kegg@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), NES, fill = Description)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("Integrative GSEA KEGG Pathways") + xlab(element_blank()) + ylab('NES') + theme(legend.position = 'none')
## padjust plot
ccle_integrate_kegg@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), NES, fill = p.adjust)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("Integrative GSEA KEGG Pathways") + xlab(element_blank()) + ylab('NES') + scale_fill_gradient(low = 'blue', high = 'red')

# Overrepresentation analysis
gene_list_filt <- gene_list[abs(gene_list) > 2]

ccle_integrate_enrich <- enricher(gene = names(gene_list_filt), universe = names(gene_list), TERM2GENE = m_t2g)
# x <- ccle_integrate_enrich@result$GeneRatio
# y <- sapply(strsplit(x, split = '/'), '[', 1)

ccle_integrate_enrich@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), Count, fill = p.adjust)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("Integrative Functional Enrichment Pathways") + xlab(element_blank()) + ylab('Gene Count') + scale_fill_gradient(low = 'blue', high = 'red')


# Reactome
ccle_integrate_reactome <- gsePathway(gene_list, organism = 'human', pvalueCutoff = 1)
ccle_integrate_reactome@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), NES, fill = Description)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("Integrative GSEA Reactome Pathways") + xlab(element_blank()) + ylab('NES') + theme(legend.position = 'none')
ccle_integrate_reactome@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), NES, fill = p.adjust)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("Integrative GSEA Reactome Pathways") + xlab(element_blank()) + ylab('NES') + scale_fill_gradient(low = 'blue', high = 'red')

'---------------------------------------------------------'
# RNA-seq DE and pathway analysis
rnaseq_data <- read.delim2(file = 'CCLE_RNAseq_genes_counts_20180929.gct', sep = '\t', skip = 2, header = T,
                           stringsAsFactors = F)
rnaseq_names <- rnaseq_data[,1:2]
line_metadata <- read.csv('breast_cell_lines_V2.csv', header = T, stringsAsFactors = F, row.names = 1)
line_metadata$ERsign <- as.factor(line_metadata$ERsign)

rnaseq_data <- rnaseq_data[,match(rownames(line_metadata), colnames(rnaseq_data))]
rownames(rnaseq_data) <- rnaseq_names$Name

all(rownames(line_metadata) == colnames(rnaseq_data))

rnaseq_DE <- DESeqDataSetFromMatrix(countData = rnaseq_data,
                                    colData = line_metadata,
                                    design = ~ERsign)
rnaseq_DE$ERsign <- relevel(rnaseq_DE$ERsign, ref = "neg")

rnaseq_DE <- DESeq(rnaseq_DE)
res <- results(rnaseq_DE, contrast = c("ERsign", "neg", "pos"), pAdjustMethod="fdr")
summary(res)

res_data <- as.data.frame(res)
res_data$gene_name <- rnaseq_names$Description[match(rownames(res_data), rnaseq_names$Name)]

# pathway analysis
rnaseq_orig_gene_list <- res_data$log2FoldChange
names(rnaseq_orig_gene_list) <- res_data$gene_name
rnaseq_orig_gene_list <- sort(rnaseq_orig_gene_list, decreasing = T)

rnaseq_orig_gene_convert <- bitr(names(rnaseq_orig_gene_list), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Hs.eg.db", drop = F)
rnaseq_orig_gene_list <- rnaseq_orig_gene_list[names(rnaseq_orig_gene_list) %in% rnaseq_orig_gene_convert$SYMBOL]
names(rnaseq_orig_gene_list) <- rnaseq_orig_gene_convert$ENTREZID 
rnaseq_orig_gene_list <- rnaseq_orig_gene_list[!is.na(rnaseq_orig_gene_convert$ENTREZID)]
hist(rnaseq_orig_gene_list)

## GSEA
rnaseq_orig_gsea <- GSEA(rnaseq_orig_gene_list, TERM2GENE = m_t2g, pvalueCutoff = .05)

## functional enrichment
rnaseq_orig_gene_list_filt <- rnaseq_orig_gene_list[abs(rnaseq_orig_gene_list) >= 2]
rnaseq_orig_enrich <- enricher(names(rnaseq_orig_gene_list_filt), universe = names(rnaseq_orig_gene_list), TERM2GENE = m_t2g)



'---------------------------------------------------------'
# RNA-seq FC only
rnaseq <- ccle_integrate %>% dplyr::select(gene_name, log2RNAFC, ENTREZID)

rnaseq_gene_list <- rnaseq$log2RNAFC
names(rnaseq_gene_list) <- rnaseq$ENTREZID
rnaseq_gene_list <- sort(rnaseq_gene_list, decreasing = T)

rnaseq_gene_list_filt <- rnaseq_gene_list[abs(rnaseq_gene_list) > 1]

# GSEA
rnaseq_gsea_C2 <- GSEA(rnaseq_gene_list, TERM2GENE = m_t2g, pvalueCutoff = 1)
rnaseq_gsea_C2@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), NES, fill = Description)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("RNA-seq GSEA Pathways") + xlab(element_blank()) + ylab('NES') + theme(legend.position = 'none')


# KEGG
rnaseq_kegg <- gseKEGG(rnaseq_gene_list, pvalueCutoff = 1)
rnaseq_kegg@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), NES, fill = Description)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("RNA-seq GSEA Pathways") + xlab(element_blank()) + ylab('NES') + theme(legend.position = 'none')


# Overrepresentation analysis
rnaseq_enrich <- enricher(names(rnaseq_gene_list_filt), universe = names(rnaseq_gene_list), TERM2GENE = m_df_filt)
rnaseq_enrich@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), GeneRatio, fill = p.adjust)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("RNA-seq Functional Enrichment Pathways") + xlab(element_blank()) + ylab('Gene Count') + scale_fill_gradient(low = 'blue', high = 'red')


# Reactome
rnaseq_reactome <- gsePathway(rnaseq_gene_list, organism = 'human', pvalueCutoff = 1)
rnaseq_reactome@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), NES, fill = Description)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("RNA-seq GSEA Pathways") + xlab(element_blank()) + ylab('NES') + theme(legend.position = 'none')


'--------------------------------------'
# Comparing pathways

'------------------------------------------------------------------------------------------------'
# Full list integrative pathway analysis
## include genes that are exclusively from RNA-seq or RPPA datasets or inclusive

# RNA-seq
rnaseq_only <- read.csv('RNA_only_CCLE_analysis_v1.csv', header = T, stringsAsFactors = F)
rnaseq_only <- rnaseq_only[,-1]
summary(rnaseq_only$log2RNAFC)
rnaseq_only <- rnaseq_only[!is.na(rnaseq_only$log2RNAFC),]
rnaseq_only <- rnaseq_only[!is.infinite(rnaseq_only$log2RNAFC),]
summary(rnaseq_only$log2RNAFC)

convert <- bitr(rnaseq_only$gene_name, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Hs.eg.db", drop = F)
rnaseq_only$ENTREZID <- convert$ENTREZID
rnaseq_only <- rnaseq_only[!is.na(rnaseq_only$ENTREZID),]

# RPPA
rppa_only <- read.csv('RPPA_only_final_CCLE_analysis_v1.csv', header = T, stringsAsFactors = F)
rppa_only <- rppa_only[-1,-1]

convert <- bitr(rppa_only$gene_name, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = "org.Hs.eg.db", drop = F)
rppa_only$ENTREZID <- convert$ENTREZID

# convert gene names to Entrez ID
ccle_all_genes <- c(ccle_integrate$max_log_FC, rnaseq_only$log2RNAF, rppa_only$log2RPPAFC)
names(ccle_all_genes) <- c(ccle_integrate$ENTREZID, rnaseq_only$ENTREZID, rppa_only$ENTREZID) 
ccle_all_genes <- sort(ccle_all_genes, decreasing = T)
summary(ccle_all_genes)

# GSEA
ccle_all_gsea <- GSEA(ccle_all_genes, TERM2GENE = m_t2g, pvalueCutoff = .05)
summary(ccle_all_gsea@result$p.adjust)
ccle_all_gsea@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), NES, fill = p.adjust)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("All-Integrative GSEA Pathways") + xlab(element_blank()) + ylab('NES') + scale_fill_gradient(low = 'blue', high = 'red')

# GSEA KEGG
ccle_all_gskegg <- gseKEGG(ccle_all_genes, pvalueCutoff = 1)
summary(ccle_all_gskegg@result$pvalue)
ccle_all_gskegg@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), NES, fill = Description)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("All-Integrative GSEA KEGG Pathways") + xlab(element_blank()) + ylab('NES') + theme(legend.position = 'none')
ccle_all_gskegg@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), NES, fill = p.adjust)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("All-Integrative GSEA KEGG Pathways") + xlab(element_blank()) + ylab('NES') + scale_fill_gradient(low = 'blue', high = 'red')

# GSEA Reactome
ccle_all_gsreactome <- gsePathway(ccle_all_genes, organism = 'human', pvalueCutoff = 1)
summary(ccle_all_gsreactome@result$pvalue)
ccle_all_gsreactome@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), NES, fill = Description)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("All-Integrative GSEA Reactome Pathways") + xlab(element_blank()) + ylab('NES') + theme(legend.position = 'none')
ccle_all_gsreactome@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), NES, fill = p.adjust)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("All-Integrative GSEA Reactome Pathways") + xlab(element_blank()) + ylab('NES') + scale_fill_gradient(low = 'blue', high = 'red')

# Functional Enrichment
ccle_all_genes_top <- ccle_all_genes[abs(ccle_all_genes) > 2]
summary(abs(ccle_all_genes_top))

ccle_all_enrich <- enricher(gene = names(ccle_all_genes_top), universe = names(ccle_all_genes), TERM2GENE = m_t2g)

ccle_all_enrich@result %>% top_n(-25, wt = pvalue) %>% ggplot(aes(reorder(Description, -pvalue), Count, fill = p.adjust)) + geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() + theme_bw() + ggtitle("All-Integrative Functional Enrichment Pathways") + xlab(element_blank()) + ylab('Gene Count') + scale_fill_gradient(low = 'blue', high = 'red')

write.csv(data.frame(gene = names(ccle_all_genes), log2FC = ccle_all_genes), file = 'ccle_all_genes.csv', row.names = F)
write.table(data.frame(gene = names(ccle_all_genes), log2FC = ccle_all_genes), file = 'ccle_all_genes.rnk', row.names = F, sep = '\t')



