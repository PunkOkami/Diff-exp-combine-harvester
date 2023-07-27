"""	DE Combine Harvester - This program is for robust differential expression analysis with single terminal command
	Copyright (C) 2023 PunkOkami

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.
	
	E-mail contact can be found on my github page: https://github.com/PunkOkami"""

import csv
from pathlib import Path
import subprocess
from matplotlib_venn import venn2
from matplotlib import pyplot as plt
import argparse
from biomart import BiomartServer
from tqdm import tqdm
import seaborn as sns
from time import sleep


def rsem_reading_data(data_dir: str) -> (dict[str: dict[str: float]], list[str], list[str]):
	# searching provided dir to find all files matching rsem output
	data_paths = Path(data_dir).rglob('*.genes.results')
	
	# reading data from files and combining them into one big object as well as reading gene names and sample names
	data_dict = {}
	gene_names = []
	for path in data_paths:
		sample_name = path.name.split('.')[0]
		in_file = open(path)
		reader = csv.reader(in_file, delimiter='\t')
		first_line = True
		name_index = 0
		counts_index = 0
		data = {}
		for line in reader:
			if first_line:
				name_index = line.index('gene_id')
				counts_index = line.index('expected_count')
				first_line = False
				continue
			gene_name = line[name_index]
			gene_count = round(float(line[counts_index]))
			data[gene_name] = gene_count
			gene_names.append(gene_name)
		data_dict[sample_name] = data
		in_file.close()
	sample_names = list(data_dict.keys())
	
	return data_dict, sample_names, gene_names


def salmon_reading_data(data_dir: str) -> (dict[str: dict[str: float]], list[str], list[str]):
	# searching provided dir to find all files matching salmon output
	data_paths = Path(data_dir).rglob('quant.genes.sf')
	
	# reading data from files and combining them into one big object as well as reading gene names and sample names
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
	sample_names = list(data_dict.keys())
	
	return data_dict, sample_names, gene_names


# argparsing
parser = argparse.ArgumentParser(prog='DE Combine Harvester',
								 description='This program is for robust differential expression analysis with single terminal command')
parser.add_argument('dirname', help='path to directory containing results of one of supported programs, seek input_type')
parser.add_argument('input_type', help='one of [salmon, rsem]. Specifies what program was used to calculate gene counts')
parser.add_argument('design_filename', help='path to file explaining experiment design, see README for information how to write one')
parser.add_argument('-o', '--out_dir', help='Optional argument used to specify what directory save results to, defaults to cwd, will be created if non-existent', default=Path.cwd(), type=Path)
parser.add_argument('-fc', '--fc_cut_off', help='Cut off for log fc values when selecting genes to be analysed by comparison, defaults to 1.5', default=1.5, type=float)
args = parser.parse_args()
data_dir = args.dirname
input_type = args.input_type
output_dir = args.out_dir
if output_dir != Path.cwd():
	output_dir.mkdir(exist_ok=True)
design_path = args.design_filename
fc_cut_off = args.fc_cut_off
scripts_wd = Path(__file__).parent

# priting copyright notice
print('DE Combine Harvester Copyright (C) 2023 PunkOkami', 'This program comes with ABSOLUTELY NO WARRANTY',
	  'This is free software, and you are welcome to redistribute it under certain conditions', 'For details see LICENSE', sep='\n')
sleep(10)

print('RUNNING DE ANALYSIS')

print('Loading data')
data_dict = {}
sample_names = []
gene_names = []
# loading data according to what program was used to count up genes
if input_type == 'rsem':
	data_dict, sample_names, gene_names = rsem_reading_data(data_dir)
elif input_type == 'salmon':
	data_dict, sample_names, gene_names = salmon_reading_data(data_dir)
else:
	print('Incorrect input_type')

# constructing sample_data file in a way for rows to be in order with colnames in data file
# also checking limiting samples to sample names in design file

# reading experiment design file, this file type is explained in detail in README
design = {}
groups = []
design_file = open(design_path)
for line in design_file:
	line = line.strip()
	line = line.split(': ')
	group = line[0]
	groups.append(group)
	samples = line[1].split(', ')
	for sample in samples:
		design[sample] = group
design_file.close()
# this variable will be used in naming output tsv files
design_name = '_vs_'.join(groups)

# pulling out list of samples in design file and checking if all of them are in sample_names
samples_in_design = list(design.keys())
for sample in samples_in_design:
	if sample not in sample_names:
		print(f'ERROR: sample {sample} file was not found in data directory, but is in experiment design file')
		exit(0)

# ordering samples dictionary according to sample_names order and creating samples_in_design list for filtering data_dict
samples_order = {sample: index for index, sample in enumerate(sample_names)}
design = sorted(design.items(), key=lambda item: samples_order[item[0]])
samples_in_design = [tup[0] for tup in design]

# creating Workdata dir just in case it does not exist yet to have some space for files needed internally
Path(scripts_wd, 'Workdata').mkdir(exist_ok=True)

# writing table saying what group is each sample in needed for R scripts
sample_data_file = open(Path(scripts_wd, 'Workdata/sample_data.tsv'), mode='w')
sample_data_writer = csv.writer(sample_data_file, delimiter='\t')
sample_data_writer.writerow(['Sample_name', 'Sample_group'])
for row in design:
	sample_data_writer.writerow(row)
sample_data_file.close()

# ordering samples to fit order in samples_in_design and filtering out ones that are not in design file
design_order = {sample_name: index for index, sample_name in enumerate(samples_in_design)}
data_dict = {sample: data for sample, data in data_dict.items() if sample in samples_in_design}
data_dict = sorted(data_dict.items(), key=lambda item: design_order[item[0]])

# putting all data in one single dictionary for easier writing to de_counts file
gene_dict = {gene: [] for gene in gene_names}
for sample, genes in data_dict:
	for gene, count in genes.items():
		gene_dict[gene].append(count)

# saving data to file to be read into R easily
de_counts_file = open(Path(scripts_wd, 'Workdata/de_counts.tsv'), mode='w')
out_writer = csv.writer(de_counts_file, delimiter='\t')
columns_names = ['Gene_ID']
columns_names.extend(samples_in_design)
out_writer.writerow(columns_names)
for gene, gene_row in gene_dict.items():
	out_row = [gene]
	out_row.extend([str(num) for num in gene_row])
	out_writer.writerow(out_row)
de_counts_file.close()

# making sure there is directory to save all Graphs
graphs_dir = Path(output_dir, 'Graphs')
graphs_dir.mkdir(exist_ok=True)
for dir_name in ['Comparison', 'DESeq2', 'edgeR']:
	Path(graphs_dir, dir_name).mkdir(exist_ok=True)

# calling R scripts and providing them all needed paths
print('Running DESeq2 analysis')
subprocess.call(f'Rscript {Path(scripts_wd, "deseq-analysis.R")} {scripts_wd} {output_dir} {fc_cut_off}', shell=True)
print('Running EdgeR analysis')
subprocess.call(f'Rscript {Path(scripts_wd, "edger-analysis.R")} {scripts_wd} {output_dir} {fc_cut_off}', shell=True)
print('DE ANALYSIS COMPLETE, LOADING RESULTS')

# loading results from DESeq2 script
deseq_results_file = open(Path(scripts_wd, 'Workdata/DESeq_results.tsv'))
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

print('\n\n')
print('Reading DESeq results')

# filtering results by fc_cut_off
deseq_results = {gene: gene_data for gene, gene_data in deseq_results.items() if abs(gene_data['log_fc']) > fc_cut_off}
deseq_results = dict(sorted(deseq_results.items(), key=lambda gene: gene[1]['fc']))
deseq_gene_ids = set(deseq_results.keys())
print(f'There is {len(deseq_results)} genes with Fold Change biologically relevant')

print('\n\n')
print('Reading EdgeR results')

# loading results from EdgeR script
edgar_results_file = open(Path(scripts_wd, 'Workdata/edger_results.tsv'))
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

# filtering results by fc_cut_off
edgar_results = {gene: gene_data for gene, gene_data in edgar_results.items() if abs(gene_data['log_fc']) > fc_cut_off}
edgar_results = dict(sorted(edgar_results.items(), key=lambda gene: gene[1]['fc']))
edgar_gene_ids = set(edgar_results.keys())
print(f'There are {len(edgar_results)} genes with Fold Change biologically relevant')

# comparing results from two methods by Venn diagram showing how many genes overlap
print('\n\n')
print('COMPARING RESULTS')
genes_in_both = list(deseq_gene_ids.intersection(edgar_gene_ids))
genes_in_both.sort()
print(f'There are {len(genes_in_both)} genes that both methods agree on')
venn2([edgar_gene_ids, deseq_gene_ids], set_labels=('EdgaR', 'DESeq2'))
plt.title('Genes found by two methods')
plt.savefig('Graphs/Comparison/venn.png', format='png')
plt.close()

# finding what biological role genes has accoriding to Biomart database

# testing whether script can connect to database, connection fails sometimes
print('Connecting to Biomart database, will take a bit, do not stop the execution')
try:
	server = BiomartServer('http://useast.ensembl.org/biomart')
except:
	print('Error occured while trying to connect to Biomart.\nHappens sometimes, try again in 15 minutes')
	exit(1)

# if connection was established, asking about genes one by one
ensembl_mart = server.datasets['hsapiens_gene_ensembl']
attributes = ['name_1006', 'namespace_1003', 'external_gene_name']
go_data = {}
gene_names = {}
print('Asking about genes')
# if there are more than 100 genes to ask about code limits traffic a bit
if len(genes_in_both) > 100:
	print('There are more than 100 genes to ask about, requests will be cut into chunks of 100 to limit traffic.')
	print('All genes will be asked about, it will just take some time')
pbar = tqdm(total=len(genes_in_both), desc='Processing genes')
for i,  gene_id in enumerate(genes_in_both):
	go_terms = []
	gene_name = ''
	response = ensembl_mart.search({'attributes': attributes, 'filters': {'ensembl_gene_id': gene_id.split('.')[0]}})
	response = response.raw.data.decode('utf-8').split('\n')[:-1]
	if len(response) != 0:
		go_terms = [line.split('\t')[0] for line in response if line.split('\t')[1] == 'biological_process'][:-1]
		gene_name = response[0].split('\t')[-1]
	if gene_name == '':
		gene_name = 'NA'
	gene_names[gene_id] = gene_name
	if len(go_terms) == 0:
		go_data[gene_id] = ['Unknown']
	else:
		go_data[gene_id] = go_terms
	pbar.update(1)
	i += 1
	if i % 100 == 0:
		print('Waiting 1 minute')
		sleep(60)
pbar.close()

# constructing catplot_data to show catplot of functions that changed the most
catplot_data = {}
for gene_id, functions in go_data.items():
	if functions == ['Unknown']:
		continue
	edgar_log_fc = edgar_results[gene_id]['log_fc']
	deseq_log_fc = deseq_results[gene_id]['log_fc']
	mean_fc = round((deseq_log_fc + edgar_log_fc) / 2, 3)
	for function in functions:
		mean_fcs = catplot_data.get(function, [])
		mean_fcs.append(mean_fc)
		catplot_data[function] = mean_fcs

# sorting data and selecting functions that have highest fold change
de_functions_file = open(Path(output_dir, f'{design_name}_de_functions.tsv'), mode='w')
de_functions_writer = csv.writer(de_functions_file, delimiter='\t')
de_functions_writer.writerow(['Biological_function', 'values_of_FC'])
catplot_data = sorted(catplot_data.items(), key=lambda tup: max(list(map(abs, tup[1]))), reverse=True)
catplot_names = []
catplot_values = []
catplot_hue = []
goal = 0

# except just filtering first 50 functions, it finds all functions that have values of fold change in first 50
# as one gene tends to be connected to many functions, this provides cut off by value and index at the same time
# also writes out tsv table with all functions
for i, tup in enumerate(catplot_data):
	function = tup[0]
	values = tup[1]
	if i <= 50 or max(list(map(abs, values))) == goal:
		catplot_values.extend(values)
		catplot_names.extend([function for i in range(len(values))])
		catplot_hue.extend([True if value > 0 else False for value in values])
	else:
		break
	if i == 50:
		goal = max(list(map(abs, values)))
	values = [str(val) for val in values]
	values = ', '.join(values)
	writer_row = [function, values]
	de_functions_writer.writerow(writer_row)
print('\n\n')
print(f'Functions result tsv table saved to {output_dir}/{design_name}_de_functions.tsv')

# plotting function expression change
catplot_data = {'names': catplot_names, 'values': catplot_values, 'Change': catplot_hue}
g = sns.catplot(catplot_data, x='values', y='names', hue='Change',  height=7, aspect=1.5, orient='v', palette='dark')
g.set_axis_labels('log FC', 'Biological function', fontsize=13, fontweight='bold')
g.figure.suptitle('Fold Change of biological functions with biggest change in expression\n\n', fontweight='bold')
g.fig.set_figwidth(11)
g.fig.set_figheight(9)
g.set(xticklabels=[])
plt.savefig('Graphs/Comparison/FC_biological_functions.png', format='png')
plt.close()

# saving data about all genes and their change to a tsv table
results_file = open(Path(output_dir, f'{design_name}_de_results.tsv'), mode='w')
results_writer = csv.writer(results_file, delimiter='\t')
genes_in_both_dict = {}
results_writer.writerow(['Gene_ID', 'Gene_name', 'FC_by_DESeq2', 'FC_by_EdgeR', 'ratio_of_FC values', 'p-adj_by_DESeq2', 'p-adj_by_EdgeR', 'biological functions'])
for gene_id in genes_in_both:
	gene_name = gene_names[gene_id]
	edgar_fc = edgar_results[gene_id]['fc']
	edgar_padj = edgar_results[gene_id]['padj']
	deseq_fc = deseq_results[gene_id]['fc']
	deseq_padj = deseq_results[gene_id]['padj']
	mean_fc = (deseq_fc+edgar_fc)/2
	ratio = deseq_fc/edgar_fc
	bio_func = go_data[gene_id]
	bio_func = ', '.join(bio_func)
	row = [gene_id, gene_name, deseq_fc, edgar_fc, ratio, deseq_padj, edgar_padj, bio_func]
	results_writer.writerow(row)
print('\n')
print(f'Genes result tsv table saved to {output_dir}/{design_name}_de_results.tsv')
