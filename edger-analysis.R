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

# perofming test to see differential expression of genes
edgar_results = exactTest(data)
edgar_results = topTags(edgar_results, n=Inf)
edgar_results = edgar_results$table
edgar_results = edgar_results[edgar_results$FDR < 0.05 & complete.cases(edgar_results$FDR), ]
edgar_results = edgar_results[order(edgar_results$FDR, decreasing = FALSE), ]
