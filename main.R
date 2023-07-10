library(tximport)

control_samples = c('SRR1747395', 'SRR1747397', 'SRR1747399')
test_samples =  c('SRR1747299', 'SRR1747301', 'SRR1747303')
samples = c(test_samples, control_samples)

files <- file.path("SALMON_OUT", list.files("SALMON_OUT"), "quant.genes.sf")
names(files) = samples
data <- tximport(files, type = "salmon", tx2gene = NULL, txIn = FALSE, geneIdCol = "Name")
