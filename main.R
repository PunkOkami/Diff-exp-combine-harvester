# loading needed libs
library(tximport)
library(DESeq2)

# setting sample names and creating factor assigning samples to groups
control_samples = c('SRR1747395', 'SRR1747397', 'SRR1747399')
test_samples =  c('SRR1747299', 'SRR1747301', 'SRR1747303')
samples = c(test_samples, control_samples)
sample_factor = gl(2, 3, labels = c('treatment', 'reference'))

# finding files and loading them into workspace
files <- file.path("SALMON_OUT", list.files("SALMON_OUT"), "quant.genes.sf")
names(files) = samples
data <- tximport(files, type = "salmon", tx2gene = NULL, txIn = FALSE, geneIdCol = "Name")

# creating colData for DESeq2
coldata = data.frame(filename = files, sample_group = sample_factor)
deseq2_data = DESeqDataSetFromTximport(data, coldata, design = ~sample_group)
dds = DESeq(deseq2_data)
resultsNames(dds)
