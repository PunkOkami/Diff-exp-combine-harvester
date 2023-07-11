# loading needed libs
library(tximport)
library(DESeq2)
library(sets)
library(pheatmap)

# setting sample names and creating factor assigning samples to groups
control_samples = c('SRR1747395', 'SRR1747397', 'SRR1747399')
test_samples =  c('SRR1747299', 'SRR1747301', 'SRR1747303')
samples = c(test_samples, control_samples)
sample_factor = gl(2, 3, labels = c('treatment', 'reference'))

# finding files and loading them into workspace
files <- file.path("SALMON_OUT", list.files("SALMON_OUT"), "quant.genes.sf")
names(files) = samples
data <- tximport(files, type = "salmon", tx2gene = NULL, txIn = FALSE, geneIdCol = "Name")

# creating colData for DESeq2 and changing object from tximport object into DESeq2 data set
coldata = data.frame(filename = files, sample_group = sample_factor)
deseq2_data = DESeqDataSetFromTximport(data, coldata, design = ~sample_group)

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
important_genes_counts = data$counts[which(rownames(data$counts) %in% small_padj_x_high_fc), ]
pheatmap(important_genes_counts, cluster_cols = FALSE)

# saving filtered results to tsv file
write.table(res, file = 'Results/DESeq_results.tsv', sep = '	', col.names = NA)