# loading needed libs
library(DESeq2)
library(sets)
library(pheatmap)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(readr)

# reading data about sample groups from file produced in Python
sample_data = read_delim('sample_data.tsv', delim = '\t', show_col_types = FALSE)
sample_data = data.frame(sample_data)
sample_names = sample_data[, 1]
sample_data = sample_data[, -1]
sample_data = data.frame(sample_group = sample_data, row.names = sample_names)
sample_data$sample_group = as.factor(sample_data$sample_group)

# loading count_data from a file produced in Python
count_data = read_delim('de_counts.tsv', delim = '\t', show_col_types = FALSE)
count_data = data.frame(count_data)
gene_ids = count_data[, 1]
count_data = count_data[, -1]
rownames(count_data) = gene_ids

# creating colData for DESeq2 and changing object from tximport object into DESeq2 data set
deseq2_data = DESeqDataSetFromMatrix(countData = round(count_data), colData = sample_data, design = ~sample_group)

# performing DESeq analysis and printing results in few ways
dds = DESeq(deseq2_data)
res = results(dds)
res =  res[complete.cases(res$padj) & res$padj < 0.05, ]
res =  res[order(res$padj), ]
summary(res)
print('Genes with smallest padj')
head(res)
print('Genes with highest absolute log2FC value')
head(res[order(abs(res$log2FoldChange), decreasing = TRUE), ])

# ploitting MA visualisation
plotMA(res, main='MA plot')

# ploting counts of 6 genes with smallest padj
small_padj = head(row.names(res), 6)
par(mfrow = c(2, 3))
for (gene in small_padj) {
  plotCounts(dds, gene = gene, intgroup = 'sample_group')
}

# looking for genes with both very high absolute fold change as well as small padj
high_fc = rownames(res[abs(res$log2FoldChange) > 1.5, ])
small_padj_x_high_fc = intersect(rownames(res), high_fc)
print('Genes that show both low p-adj and high FC')
small_padj_x_high_fc

# plotting vulcanoplot showing relation between FC and padj and pointing out genes found using sets
par(mfrow = c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch = 20, main = "Volcano plot", col = 'blue'))
with(subset(res, padj < 0.05 & abs(log2FoldChange) > 1.5), points(log2FoldChange, -log10(pvalue), pch = 20, col = "red"))

# creating heatmap of counts of genes with high FC
important_genes_counts = count_data[which(rownames(count_data) %in% small_padj_x_high_fc), ]
pheatmap(important_genes_counts, cluster_cols = FALSE, color = hcl.colors(50, "plasma"))

# performing and ploting PCA
rld = rlog(dds, blind = FALSE)
plotPCA(rld, intgroup = 'sample_group') + geom_text(aes(label = sample_names), vjust = 2, size = 3)

# calculating sample-to-sample distances
distsRL = dist(t(assay(rld)))
dist_mat = as.matrix(distsRL)
hc = hclust(distsRL)
hmcol = colorRampPalette(brewer.pal(9, 'PRGn'))(100)
heatmap.2(dist_mat, Rowv = as.dendrogram(hc), symm = TRUE, trace = "none", col = rev(hmcol), margin = c(13, 13))

# saving filtered results to tsv file
write.table(res, file = 'Results/DESeq_results.tsv', sep = '\t', col.names = NA)