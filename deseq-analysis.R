    # DE Combine Harvester - This program is for robust differential expression analysis with single terminal command
    # Copyright (C) 2023 PunkOkami
    
    # This program is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 3 of the License, or
    # (at your option) any later version.
    
    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.
    
    # You should have received a copy of the GNU General Public License
    # along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
    # E-mail contact can be found on my github page: https://github.com/PunkOkami"

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(sets))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(readr))

# parsing args
args = commandArgs(trailingOnly = TRUE)
workspace = args[1]
output_dir = args[2]
# this arg is cut off set by Python input
fc_cut_off = as.numeric(args[3])

# setting wd to where program is located
setwd(workspace)

# reading data about sample groups from file produced in Python
sample_data = suppressMessages(read_delim('Workdata/sample_data.tsv', delim = '\t', show_col_types = FALSE))
sample_data = data.frame(sample_data)
sample_names = sample_data[, 1]
sample_data = sample_data[, -1]
sample_data = data.frame(sample_group = sample_data, row.names = sample_names)
sample_data$sample_group = as.factor(sample_data$sample_group)

# loading count_data from a file produced in Python
count_data = suppressMessages(read_delim('Workdata/de_counts.tsv', delim = '\t', show_col_types = FALSE))
count_data = data.frame(count_data)
gene_ids = count_data[, 1]
count_data = count_data[, -1]
rownames(count_data) = gene_ids

# using dataset read from a file and sample data to create DESeqDataSet
deseq2_data = suppressMessages(DESeqDataSetFromMatrix(countData = round(count_data), colData = sample_data, design = ~sample_group))

# performing DESeq analysis and filtering results by p-adj
dds = DESeq(deseq2_data, quiet = TRUE)
res = results(dds)
res =  res[complete.cases(res$padj) & res$padj < 0.05, ]
res =  res[order(res$padj), ]

# settiong wd to output_dir to put graphs where needed
setwd(output_dir)

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
high_fc = rownames(res[abs(res$log2FoldChange) > fc_cut_off, ])
small_padj_x_high_fc = intersect(rownames(res), high_fc)

# plotting vulcanoplot showing relation between FC and padj
png(filename = 'Graphs/DESeq2/Vulcano.png')
with(res, plot(log2FoldChange, -log10(pvalue), pch = 20, main = "Volcano plot", col = 'blue'))
with(subset(res, padj < 0.05 & abs(log2FoldChange) > fc_cut_off), points(log2FoldChange, -log10(pvalue), pch = 20, col = "red"))
legend(-8, 220,  legend = c('absolute FC < fc_cut_off', 'absolute FC > fc_cut_off'), col = c('blue', 'red'))
garbage = dev.off()

# creating heatmap of counts of genes with high FC
important_genes_counts = count_data[which(rownames(count_data) %in% small_padj_x_high_fc), ]
important_genes_counts = sweep(important_genes_counts, 2, sizeFactors(dds), '/')
if (nrow(important_genes_counts) <= 75) {
  fontsize = round(10/(nrow(important_genes_counts)%%15), digits = 2)
  if (fontsize == 0){fontsize = 10}
  pheatmap(important_genes_counts, cluster_cols = FALSE, fontsize_row = fontsize,
         color = hcl.colors(50, "plasma"), filename = 'Graphs/DESeq2/heatmap.png')
} else {print("Heatmap of counts for genes with high fc would be unreadable, won't be plotted")}


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

# jumping back to where code is placed to save workdata properly
setwd(workspace)

res = cbind(data.frame(GeneID = rownames(res), Fold_change=2^res$log2FoldChange), data.frame(res))
# saving filtered results to tsv file
write.table(res, file = 'Workdata/DESeq_results.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
