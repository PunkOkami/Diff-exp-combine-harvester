library(edgeR)
library(tximport)

# setting sample names and creating factor assigning samples to groups
control_samples = c('SRR1747395', 'SRR1747397', 'SRR1747399')
test_samples =  c('SRR1747299', 'SRR1747301', 'SRR1747303')
samples = c(test_samples, control_samples)
sample_factor = gl(2, 3, labels = c('treatment', 'reference'))

# finding files and loading them into workspace
files <- file.path("SALMON_OUT", list.files("SALMON_OUT"), "quant.genes.sf")
names(files) = samples
data <- tximport(files, type = "salmon", tx2gene = NULL, txIn = FALSE, geneIdCol = "Name")

# transforms data into format used by edgeR and performs basic
data = DGEList(data$counts, group = sample_factor)
data = calcNormFactors(data)
data = estimateCommonDisp(data)
data = estimateTagwiseDisp(data)

# plots MDS plot to show how samples compare
plotMDS(data, method="bcv", col=as.numeric(data$samples$group))
legend("bottomleft", as.character(unique(data$samples$group)), col=1:2, pch=20)

# perofming test to see differential expression of genes
edgar_results = exactTest(data)

# showing summary of differential expression
print('Summary of gene expression')
de_genes = decideTestsDGE(edgar_results, p.value = 0.05, lfc = 1.5)
summary(de_genes)

# geting results of de
edgar_results = topTags(edgar_results, n=Inf)
edgar_results = edgar_results$table
edgar_results = edgar_results[edgar_results$FDR < 0.05 & complete.cases(edgar_results$FDR), ]
edgar_results = edgar_results[order(edgar_results$FDR, decreasing = FALSE), ]

# saving data to a file
edgar_results =  cbind(data.frame(Geneid=rownames(edgar_results), Fold_change=2**edgar_results$logFC), edgar_results)
write.table(edgar_results, file = "Results/edger_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
