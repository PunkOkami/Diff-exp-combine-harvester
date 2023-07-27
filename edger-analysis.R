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

# loading needed libs
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(pheatmap))
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

# transforms data into format used by edgeR and performs basic calculations needed for later analysis
edgar_data = DGEList(count_data, group = sample_data$sample_group)
edgar_data = calcNormFactors(edgar_data)
edgar_data = estimateCommonDisp(edgar_data)
edgar_data = estimateTagwiseDisp(edgar_data)

# settiong wd to output_dir to put graphs where needed
setwd(output_dir)

# plots MDS plot to show how samples compare
png(filename = 'Graphs/edgeR/MDS_plot.png')
plotMDS(edgar_data, col=as.numeric(edgar_data$samples$group), main = 'Sample similarity')
legend("bottomleft", as.character(unique(edgar_data$samples$group)), col=1:2, pch=20)
garbage = dev.off()

# perofming test to see differential expression of genes
edgar_results = exactTest(edgar_data)

# showing summary of differential expression
de_genes = decideTestsDGE(edgar_results, lfc = fc_cut_off)
important_genes = rownames(de_genes)[which(de_genes != 0)]
png(filename = 'Graphs/edgeR/Mean_difference_plot.png')
plotSmear(edgar_results, de.tags = important_genes, main = 'Mean difference plot')
garbage = dev.off()

# creating heatmap of counts of genes with high FC
important_genes_counts = count_data[which(rownames(count_data) %in% important_genes), ]
if (nrow(important_genes_counts) <= 75) {
  fontsize = round(10/(nrow(important_genes_counts)%%15), digits = 2)
  if (fontsize == 0){fontsize = 10}
  pheatmap(important_genes_counts, cluster_cols = FALSE, fontsize_row = fontsize,
         color = hcl.colors(50, "plasma"), filename = 'Graphs/edgeR/counts_heatmap.png')
} else {print("Heatmap of counts for genes with high fc would be unreadable, won't be plotted")}

# geting results of de
edgar_results = topTags(edgar_results, n=Inf)
edgar_results = edgar_results$table
edgar_results = edgar_results[edgar_results$FDR < 0.05 & complete.cases(edgar_results$FDR), ]
edgar_results = edgar_results[order(edgar_results$FDR, decreasing = FALSE), ]

# jumping back to where code is placed to save workdata properly
setwd(workspace)

# saving data to a file
edgar_results =  cbind(data.frame(GeneID=rownames(edgar_results), Fold_change=2^edgar_results$logFC), edgar_results)
write.table(edgar_results, file = "Workdata/edger_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
