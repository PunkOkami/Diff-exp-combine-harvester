import csv
from pathlib import Path
from pympler.asizeof import asizeof


def salmon_reading_data(data_dir: str):
	# Finding data files
	data_paths = Path(data_dir).rglob('quant.genes.sf')
	
	# Reading data from files and combining them into one big object
	data_dict = {}
	gene_names = []
	for path in data_paths:
		sample_name = path.parts[-2]
		reader = csv.reader(open(path), delimiter='\t')
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
	
	# Putting all data in one single matrix
	sample_names = []
	gene_dict = {gene: [] for gene in gene_names}
	
	for sample, genes in data_dict.items():
		sample_names.append(sample)
		for gene, count in genes.items():
			gene_dict[gene].append(count)
	
	out_writer = csv.writer(open('de_counts', mode='w'), delimiter='\t')
	columns_names = ['Gene ID']
	columns_names.extend(sample_names)
	out_writer.writerow(columns_names)
	

# Directory of data, later to be changed for agrpass
salmon_data_dir = 'Example_data'

salmon_reading_data(salmon_data_dir)
