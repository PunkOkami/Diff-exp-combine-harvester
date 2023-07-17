import csv
from pathlib import Path
import subprocess
from matplotlib_venn import venn2
from matplotlib import pyplot as plt


def salmon_reading_data(data_dir: str):
	# Finding data files
	data_paths = Path(data_dir).rglob('quant.genes.sf')
	
	# Reading data from files and combining them into one big object
	data_dict = {}
	gene_names = []
	for path in data_paths:
		sample_name = path.parts[-2]
		in_file = open(path)
		reader = csv.reader(in_file, delimiter='\t')
		first_line = True
		name_index = 0
		counts_index = 0
		data = {}
		for line in reader:
			if first_line:
				name_index = line.index('Name')
				counts_index = line.index('NumReads')
				first_line = False
				continue
			gene_name = line[name_index]
			gene_count = float(line[counts_index])
			data[gene_name] = gene_count
			gene_names.append(gene_name)
		data_dict[sample_name] = data
		in_file.close()
	
	# Putting all data in one single dictionary
	sample_names = []
	gene_dict = {gene: [] for gene in gene_names}
	for sample, genes in data_dict.items():
		sample_names.append(sample)
		for gene, count in genes.items():
			gene_dict[gene].append(count)
	
	# Saving data to file to be read into R easily
	output_file = open('de_counts.tsv', mode='w')
	out_writer = csv.writer(output_file, delimiter='\t')
	columns_names = ['Gene_ID']
	columns_names.extend(sample_names)
	out_writer.writerow(columns_names)
	for gene, gene_row in gene_dict.items():
		out_row = [gene]
		out_row.extend([str(num) for num in gene_row])
		out_writer.writerow(out_row)
	output_file.close()
	
	# Creating table saying what group is what sample is in
	sample_groups = ['treatment', 'treatment', 'reference', 'reference', 'reference', 'treatment']
	sample_data = []
	for sample, group in zip(sample_names, sample_groups):
		sample_data.append([sample, group])
	sample_data_file = open('sample_data.tsv', mode='w')
	sample_data_writer = csv.writer(sample_data_file, delimiter='\t')
	sample_data_writer.writerow(['Sample_name', 'Sample_group'])
	for row in sample_data:
		sample_data_writer.writerow(row)
	sample_data_file.close()


# Directory of data, later to be changed for agrpass
salmon_data_dir = 'Example_data'
print('Loading data')
# salmon_reading_data(salmon_data_dir)

# calling R scripts
print('Running DESeq2 analysis')
# subprocess.call('Rscript deseq-analysis.R', shell=True)
print('Running EdgeR analysis')
# subprocess.call('Rscript edger-analysis.R', shell=True)
print('DE analysis complete, loading results')

# loading results from DESeq2 script and doing subtle result analysis
deseq_results_file = open('Results/DESeq_results.tsv')
deseq_reader = csv.reader(deseq_results_file, delimiter='\t')
first_line = True
deseq_results = []
gene_id_index = 0
fc_index = 0
padj_index = 0
for line in deseq_reader:
	if first_line:
		gene_id_index = line.index('GeneID')
		fc_index = line.index('log2FoldChange')
		padj_index = line.index('padj')
		first_line = False
		continue
	gene_id = line[gene_id_index]
	line = [float(line[ind]) for ind in range(len(line)) if ind != gene_id_index]
	line.insert(0, gene_id)
	line = tuple(line)
	deseq_results.append(line)
deseq_results_file.close()

print('Analysing DESeq results')
print('\n\n')
print('6 genes with smallest p-adj:')
print(f'Gene ID - p-adj')
for i, gene in enumerate(deseq_results):
	if i > 5:
		break
	print(f'{gene[0]} - {gene[padj_index]}')
print('\n\n')

deseq_results = [gene for gene in deseq_results if abs(gene[fc_index]) > 1.5]
deseq_results = sorted(deseq_results, key=lambda gene: gene[fc_index])
print(f'There is {len(deseq_results)} genes with Fold Change biologically relevant')
print('Gene ID - Fold Change - p-adj')
deseq_gene_ids = set()
for gene in deseq_results:
	deseq_gene_ids.add(gene[0])
	print(f'{gene[0]} - {2**gene[fc_index]} - {gene[padj_index]}')
print('\n\n')

# loading results from EdgeR script and doing subtle result analysis
edgar_results_file = open('Results/edger_results.tsv')
edgar_reader = csv.reader(edgar_results_file, delimiter='\t')
first_line = True
edgar_results = []
gene_id_index = 0
fc_index = 0
padj_index = 0
log_fc_index = 0
for line in edgar_reader:
	if first_line:
		gene_id_index = line.index('GeneID')
		fc_index = line.index('Fold_change')
		padj_index = line.index('FDR')
		log_fc_index = line.index('logFC')
		first_line = False
		continue
	gene_id = line[gene_id_index]
	line = [float(line[ind]) for ind in range(len(line)) if ind != gene_id_index]
	line.insert(0, gene_id)
	line = tuple(line)
	edgar_results.append(line)
edgar_results_file.close()

print('Analysing EdgeR results')
print('\n\n')
print('6 genes with smallest p-adj:')
print(f'Gene ID - p-adj')
for i, gene in enumerate(edgar_results):
	if i > 5:
		break
	print(f'{gene[0]} - {gene[padj_index]}')
print('\n\n')

edgar_results = [gene for gene in edgar_results if abs(gene[log_fc_index]) > 1.5]
edgar_results = sorted(edgar_results, key=lambda gene: gene[fc_index])
print(f'There is {len(edgar_results)} genes with Fold Change biologically relevant')
print('Gene ID - Fold Change - p-adj')
edgar_gene_ids = set()
for gene in edgar_results:
	edgar_gene_ids.add(gene[0])
	print(f'{gene[0]} - {gene[fc_index]} - {gene[padj_index]}')
print('\n\n')

# comparing results from two methods
print('Comapring results')
genes_in_both = deseq_gene_ids.intersection(edgar_gene_ids)
print(len(genes_in_both))
only_deseq = deseq_gene_ids.difference(edgar_gene_ids)
only_edgar = edgar_gene_ids.difference(deseq_gene_ids)
venn2(subsets=(len(only_edgar), len(only_deseq), len(genes_in_both)), set_labels=('EdgaR', 'DESeq2'))
plt.title('Genes found by two methods')
plt.show()

print('Genes two methods agree on:')
print(' ,'.join(genes_in_both))
