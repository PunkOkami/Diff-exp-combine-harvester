# loading needed libs
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(readr))

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

# transforms data into format used by edgeR and performs basic calculations needed for later analysis
edgar_data = DGEList(count_data, group = sample_data$sample_group)
edgar_data = calcNormFactors(edgar_data)
edgar_data = estimateCommonDisp(edgar_data)
edgar_data = estimateTagwiseDisp(edgar_data)

# plots MDS plot to show how samples compare
png(filename = 'Graphs/edgeR/MDS_plot.png')
plotMDS(edgar_data, col=as.numeric(edgar_data$samples$group), main = 'Sample similarity')
legend("bottomleft", as.character(unique(edgar_data$samples$group)), col=1:2, pch=20)
garbage = dev.off()

# perofming test to see differential expression of genes
edgar_results = exactTest(edgar_data)

# showing summary of differential expression
de_genes = decideTestsDGE(edgar_results, lfc = 1.5)
important_genes = rownames(de_genes)[which(de_genes != 0)]
png(filename = 'Graphs/edgeR/Mean_difference_plot.png')
plotSmear(edgar_results, de.tags = important_genes, main = 'Mean difference plot')
garbage = dev.off()

important_genes_counts = count_data[which(rownames(count_data) %in% important_genes), ]
pheatmap(important_genes_counts, cluster_cols = FALSE,
		 color = hcl.colors(50, "cyan-magenta"), filename = 'Graphs/edgeR/counts_heatmap.png')

# geting results of de
edgar_results = topTags(edgar_results, n=Inf)
edgar_results = edgar_results$table
edgar_results = edgar_results[edgar_results$FDR < 0.05 & complete.cases(edgar_results$FDR), ]
edgar_results = edgar_results[order(edgar_results$FDR, decreasing = FALSE), ]

# saving data to a file
edgar_results =  cbind(data.frame(GeneID=rownames(edgar_results), Fold_change=2^edgar_results$logFC), edgar_results)
write.table(edgar_results, file = "Workdata/edger_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
