import csv
from pathlib import Path
import subprocess
from matplotlib_venn import venn2
from matplotlib import pyplot as plt
import argparse


def salmon_reading_data(data_dir: str):
	print('Loading data')
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
	output_file = open('Workdata/de_counts.tsv', mode='w')
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
	sample_data_file = open('Workdata/sample_data.tsv', mode='w')
	sample_data_writer = csv.writer(sample_data_file, delimiter='\t')
	sample_data_writer.writerow(['Sample_name', 'Sample_group'])
	for row in sample_data:
		sample_data_writer.writerow(row)
	sample_data_file.close()


# argparsing
parser = argparse.ArgumentParser(prog='Differential Expression Combine Harvester',
								description='This program is for robust differential expression analysis with one terminal command')
parser.add_argument('dirname', help='path to directory containing results of one of supported programs, seek input_type')
parser.add_argument('input_type', help='one of [salmon, rsem]. Specifies what program was used to calculate gene counts')
parser.add_argument('-o', '--out_file', help='Optional argument used to specify where save results in tsv', default='de_results.tsv')
args = parser.parse_args()
data_dir = args.dirname
input_type = args.input_type
output_file = args.out_file

# loading data
if input_type == 'rsem':
	print('This feature is not supported yet')
	exit(0)
elif input_type == 'salmon':
	salmon_reading_data(data_dir)
else:
	print('Incorrect input_type')

# calling R scripts
print('Running DESeq2 analysis')
subprocess.call('Rscript deseq-analysis.R', shell=True)
print('Running EdgeR analysis')
subprocess.call('Rscript edger-analysis.R', shell=True)
print('DE analysis complete, loading results')

# loading results from DESeq2 script and doing subtle result analysis
deseq_results_file = open('Workdata/DESeq_results.tsv')
deseq_reader = csv.reader(deseq_results_file, delimiter='\t')
first_line = True
deseq_results = {}
gene_id_index = 0
fc_index = 0
padj_index = 0
log_fc_index = 0
for line in deseq_reader:
	if first_line:
		gene_id_index = line.index('GeneID')
		log_fc_index = line.index('log2FoldChange')
		fc_index = line.index('Fold_change')
		padj_index = line.index('padj')
		first_line = False
		continue
	gene_id = line[gene_id_index]
	gene_data = {'padj': float(line[padj_index]), 'fc': float(line[fc_index]), 'log_fc': float(line[log_fc_index])}
	deseq_results[gene_id] = gene_data
deseq_results_file.close()

print('Analysing DESeq results')
print('\n\n')
print('6 genes with smallest p-adj:')
print(f'Gene ID - p-adj')
for i, (gene_id, gene_data) in enumerate(deseq_results.items()):
	if i > 5:
		break
	print(f'{gene_id} - {gene_data["padj"]}')
print('\n\n')

deseq_results = {gene: gene_data for gene, gene_data in deseq_results.items() if abs(gene_data['log_fc']) > 1.5}
deseq_results = dict(sorted(deseq_results.items(), key=lambda gene: gene[1]['fc']))
print(f'There is {len(deseq_results)} genes with Fold Change biologically relevant')
print('Gene ID - Fold Change - p-adj')
deseq_gene_ids = set()
for gene, gene_data in deseq_results.items():
	deseq_gene_ids.add(gene)
	print(f'{gene} - {gene_data["fc"]} - {gene_data["padj"]}')
print('\n\n')

# loading results from EdgeR script and doing subtle result analysis
edgar_results_file = open('Workdata/edger_results.tsv')
edgar_reader = csv.reader(edgar_results_file, delimiter='\t')
first_line = True
edgar_results = {}
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
	gene_data = {'padj': float(line[padj_index]), 'fc': float(line[fc_index]), 'log_fc': float(line[log_fc_index])}
	edgar_results[gene_id] = gene_data
edgar_results_file.close()

print('Analysing EdgeR results')
print('\n\n')
print('6 genes with smallest p-adj:')
print(f'Gene ID - p-adj')
for i, (gene, gene_data) in enumerate(edgar_results.items()):
	if i > 5:
		break
	print(f'{gene} - {gene_data["padj"]}')
print('\n\n')

edgar_results = {gene: gene_data for gene, gene_data in edgar_results.items() if abs(gene_data['log_fc']) > 1.5}
edgar_results = dict(sorted(edgar_results.items(), key=lambda gene: gene[1]['fc']))
print(f'There are {len(edgar_results)} genes with Fold Change biologically relevant')
print('Gene ID - Fold Change - p-adj')
edgar_gene_ids = set()
for gene, gene_data in edgar_results.items():
	edgar_gene_ids.add(gene)
	print(f'{gene} - {gene_data["fc"]} - {gene_data["padj"]}')
print('\n\n')

# comparing results from two methods by Venn diagram showing how many genes overlap
print('Comapring results')
genes_in_both = deseq_gene_ids.intersection(edgar_gene_ids)
venn2([edgar_gene_ids, deseq_gene_ids], set_labels=('EdgaR', 'DESeq2'))
plt.title('Genes found by two methods')
plt.savefig('Graphs/Comparison/venn.png', format='png')
plt.close()

# priting info about fc of genes two methods agree on and saving that data to a file
results_file = open(output_file, mode='w')
results_writer = csv.writer(results_file, delimiter='\t')
genes_in_both_dict = {}
results_writer.writerow(['GeneID', 'FC_by_DESeq2', 'FC_by_EdgeR', 'ratio_of_FC values', 'p-adj_by_DESeq2', 'p-adj_by_EdgeR'])
print('Genes both methods agree on:')
print('Gene ID - FC by DESeq2 - FC by EdgeR - ratio between DESeq2 and EdgeR values')
for gene_id in genes_in_both:
	edgar_fc = edgar_results[gene_id]['fc']
	edgar_padj = edgar_results[gene_id]['padj']
	deseq_fc = deseq_results[gene_id]['fc']
	deseq_padj = deseq_results[gene_id]['padj']
	mean_fc = (deseq_fc+edgar_fc)/2
	ratio = deseq_fc/edgar_fc
	print(f'{gene_id} - {deseq_fc} - {edgar_fc} - {ratio}')
	row = [gene_id, deseq_fc, edgar_fc, ratio, deseq_padj, edgar_padj]
	results_writer.writerow(row)
print('\n\n')
print(f'Result tsv table saved to {output_file}')
