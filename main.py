import csv
from pathlib import Path
import subprocess
from matplotlib_venn import venn2
from matplotlib import pyplot as plt
import argparse
from biomart import BiomartServer
from tqdm import tqdm
import seaborn as sns


def salmon_reading_data(data_dir: str) -> (dict[str: dict[str: float]], list[str], list[str]):
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
	sample_names = list(data_dict.keys())
	return data_dict, sample_names, gene_names


# argparsing
parser = argparse.ArgumentParser(prog='Differential Expression Combine Harvester',
								description='This program is for robust differential expression analysis with one terminal command')
parser.add_argument('dirname', help='path to directory containing results of one of supported programs, seek input_type')
parser.add_argument('input_type', help='one of [salmon, rsem]. Specifies what program was used to calculate gene counts')
parser.add_argument('design_filename', help='path to file explaining experiment design, see README for information how to write one')
parser.add_argument('-o', '--out_file', help='Optional argument used to specify where save results in tsv', default='de_results.tsv')
args = parser.parse_args()
data_dir = args.dirname
input_type = args.input_type
output_file = args.out_file
design_path = args.design_filename

# Creating table saying what group is what sample is in
print('Loading data')
data_dict = {}
sample_names = []
gene_names = []
# loading data
if input_type == 'rsem':
	print('This feature is not supported yet')
	exit(0)
elif input_type == 'salmon':
	data_dict, sample_names, gene_names = salmon_reading_data(data_dir)
else:
	print('Incorrect input_type')
	
# constructing sample_data file in a way for rows to be in order with colnames in data file
# also checking limiting samples to sample names in design file
# reading experiment design file
design = {}
design_file = open(design_path)
for line in design_file:
	line = line.strip()
	line = line.split(': ')
	group = line[0]
	samples = line[1].split(', ')
	for sample in samples:
		design[sample] = group
design_file.close()

# pulling out list of samples in design file and checking if all of them are in sample_names
samples_in_design = list(design.keys())
for sample in samples_in_design:
	if sample not in sample_names:
		print(f'ERROR: sample {sample} file was not found in data directory, but is in experiment design file')
		exit(0)

# ordering dictionary according to sample_names
samples_order = {sample: index for index, sample in enumerate(sample_names)}
design = sorted(design.items(), key=lambda item: samples_order[item[0]])
samples_in_design = [tup[0] for tup in design]

# creating table saying what group is what sample is in
sample_data_file = open('Workdata/sample_data.tsv', mode='w')
sample_data_writer = csv.writer(sample_data_file, delimiter='\t')
sample_data_writer.writerow(['Sample_name', 'Sample_group'])
for row in design:
	sample_data_writer.writerow(row)
sample_data_file.close()

# Ordering samples to fit order in samples_in_design
design_order = {sample_name: index for index, sample_name in enumerate(samples_in_design)}
data_dict = {sample: data for sample, data in data_dict.items() if sample in samples_in_design}
data_dict = sorted(data_dict.items(), key=lambda item: design_order[item[0]])

# Putting all data in one single dictionary
gene_dict = {gene: [] for gene in gene_names}
for sample, genes in data_dict:
	for gene, count in genes.items():
		gene_dict[gene].append(count)

# Saving data to file to be read into R easily
de_counts_file = open('Workdata/de_counts.tsv', mode='w')
out_writer = csv.writer(de_counts_file, delimiter='\t')
columns_names = ['Gene_ID']
columns_names.extend(sample_names)
out_writer.writerow(columns_names)
for gene, gene_row in gene_dict.items():
	out_row = [gene]
	out_row.extend([str(num) for num in gene_row])
	out_writer.writerow(out_row)
de_counts_file.close()

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
print('COMPARING RESULTS')
genes_in_both = list(deseq_gene_ids.intersection(edgar_gene_ids))
genes_in_both.sort()
venn2([edgar_gene_ids, deseq_gene_ids], set_labels=('EdgaR', 'DESeq2'))
plt.title('Genes found by two methods')
plt.savefig('Graphs/Comparison/venn.png', format='png')
plt.close()

# finding what biological role genes has accoriding to Biomart database
print('Connecting to Biomart database, will take a bit, do not stop the execution')
try:
	server = BiomartServer('http://useast.ensembl.org/biomart')
except:
	print('Error occured while trying to connect to Biomart.\nHappens sometimes, try again in 15 minutes')
	exit(1)
ensembl_mart = server.datasets['hsapiens_gene_ensembl']
attributes = ['name_1006', 'namespace_1003', 'external_gene_name']
go_data = {}
gene_names = {}
print('Asking about genes')
pbar = tqdm(total=len(genes_in_both), desc='Processing genes')
for gene_id in genes_in_both:
	response = ensembl_mart.search({'attributes': attributes, 'filters': {'ensembl_gene_id': gene_id}})
	response = response.raw.data.decode('utf-8').split('\n')[:-1]
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
pbar.close()

# constructing barplot_data
barplot_names = []
barplot_values = []
for gene_id, functions in go_data.items():
	if functions == ['Unknown']:
		continue
	edgar_fc = edgar_results[gene_id]['fc']
	deseq_fc = deseq_results[gene_id]['fc']
	mean_fc = round((deseq_fc + edgar_fc) / 2, 3)
	for function in functions:
		gene_name = gene_names[gene_id]
		barplot_names.append(function)
		barplot_values.append(mean_fc)

g = sns.catplot(x=barplot_values, y=barplot_names, height=7, aspect=1.5, orient='v')
g.set_axis_labels('FC', 'Biological function', fontsize=13, fontweight='bold')
g.figure.suptitle('Fold Change of each biological function that changed expression\n\n', fontweight='bold')
g.fig.set_figwidth(11)
g.fig.set_figheight(9)
plt.savefig('Graphs/Comparison/FC_biological_functions.png', format='png')
plt.close()

# priting info about fc of genes two methods agree on and saving that data to a file
results_file = open(output_file, mode='w')
results_writer = csv.writer(results_file, delimiter='\t')
genes_in_both_dict = {}
results_writer.writerow(['Gene_ID', 'Gene_name', 'FC_by_DESeq2', 'FC_by_EdgeR', 'ratio_of_FC values', 'p-adj_by_DESeq2', 'p-adj_by_EdgeR', 'biological functions'])
print('Genes both methods agree on:')
print('Gene ID - FC by DESeq2 - FC by EdgeR - ratio between DESeq2 and EdgeR values')
for gene_id in genes_in_both:
	gene_name = gene_names[gene_id]
	edgar_fc = edgar_results[gene_id]['fc']
	edgar_padj = edgar_results[gene_id]['padj']
	deseq_fc = deseq_results[gene_id]['fc']
	deseq_padj = deseq_results[gene_id]['padj']
	mean_fc = (deseq_fc+edgar_fc)/2
	ratio = deseq_fc/edgar_fc
	print(f'{gene_id} - {deseq_fc} - {edgar_fc} - {ratio}')
	bio_func = go_data[gene_id]
	bio_func = ', '.join(bio_func)
	row = [gene_id, gene_name, deseq_fc, edgar_fc, ratio, deseq_padj, edgar_padj, bio_func]
	results_writer.writerow(row)
print('\n\n')
print(f'Result tsv table saved to {output_file}')
