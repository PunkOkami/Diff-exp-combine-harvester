# loading needed libs
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(sets))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(readr))

# reading data about sample groups from file produced in Python
sample_data = suppressMessages(read_delim('sample_data.tsv', delim = '\t', show_col_types = FALSE))
sample_data = data.frame(sample_data)
sample_names = sample_data[, 1]
sample_data = sample_data[, -1]
sample_data = data.frame(sample_group = sample_data, row.names = sample_names)
sample_data$sample_group = as.factor(sample_data$sample_group)

# loading count_data from a file produced in Python
count_data = suppressMessages(read_delim('de_counts.tsv', delim = '\t', show_col_types = FALSE))
count_data = data.frame(count_data)
gene_ids = count_data[, 1]
count_data = count_data[, -1]
rownames(count_data) = gene_ids

# creating colData for DESeq2 and changing object from tximport object into DESeq2 data set
deseq2_data = suppressMessages(DESeqDataSetFromMatrix(countData = round(count_data), colData = sample_data, design = ~sample_group))

# performing DESeq analysis and printing results in few ways
dds = DESeq(deseq2_data, quiet = TRUE)
res = results(dds)
res =  res[complete.cases(res$padj) & res$padj < 0.05, ]
res =  res[order(res$padj), ]

# ploitting MA visualisation
png(filename = 'Graphs/DESeq2/MA_plot.png')
plotMA(res, main='MA plot')
garbage = dev.off()

# ploting counts of 6 genes with smallest padj
small_padj = head(row.names(res), 6)
png(filename = 'Graphs/DESeq2/plotCounts.png', width = 600)
par(mfrow = c(2, 3))
for (gene in small_padj) {
  plotCounts(dds, gene = gene, intgroup = 'sample_group')
}
garbage = dev.off()

# looking for genes with both very high absolute fold change as well as small padj
high_fc = rownames(res[abs(res$log2FoldChange) > 1.5, ])
small_padj_x_high_fc = intersect(rownames(res), high_fc)

# plotting vulcanoplot showing relation between FC and padj and pointing out genes found using sets
png(filename = 'Graphs/DESeq2/Vulcano.png')
with(res, plot(log2FoldChange, -log10(pvalue), pch = 20, main = "Volcano plot", col = 'blue'))
with(subset(res, padj < 0.05 & abs(log2FoldChange) > 1.5), points(log2FoldChange, -log10(pvalue), pch = 20, col = "red"))
garbage = dev.off()

# creating heatmap of counts of genes with high FC
important_genes_counts = count_data[which(rownames(count_data) %in% small_padj_x_high_fc), ]
important_genes_counts = sweep(important_genes_counts, 2, sizeFactors(dds), '/')
pheatmap(important_genes_counts, cluster_cols = FALSE,
         color = hcl.colors(50, "plasma"), filename = 'Graphs/DESeq2/heatmap.png')

# performing and ploting PCA
png(filename = 'Graphs/DESeq2/PCA_plot.png')
rld = rlog(dds, blind = FALSE)
plotPCA(rld, intgroup = 'sample_group') + geom_text(aes(label = sample_names), vjust = 2, size = 3)
garbage = dev.off()

# calculating sample-to-sample distances
distsRL = dist(t(assay(rld)))
dist_mat = as.matrix(distsRL)
hc = hclust(distsRL)
hmcol = colorRampPalette(brewer.pal(9, 'PRGn'))(100)
png(filename = 'Graphs/DESeq2/sample_to_sample_heatmap.png')
heatmap.2(dist_mat, Rowv = as.dendrogram(hc), symm = TRUE, trace = "none", col = rev(hmcol), margin = c(13, 13))
garbage = dev.off()

res = cbind(data.frame(GeneID = rownames(res)), data.frame(res))
# saving filtered results to tsv file
write.table(res, file = 'Results/DESeq_results.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
