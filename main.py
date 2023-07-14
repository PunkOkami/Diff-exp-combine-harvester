import csv
from pathlib import Path


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
	output_file = open('Example_data/de_counts.tsv', mode='w')
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
	sample_data_file = open('Example_data/sample_data.tsv', mode='w')
	sample_data_writer = csv.writer(sample_data_file, delimiter='\t')
	sample_data_writer.writerow(['Sample_name', 'Sample_group'])
	for row in sample_data:
		sample_data_writer.writerow(row)
	sample_data_file.close()


# Directory of data, later to be changed for agrpass
salmon_data_dir = 'Example_data'
salmon_reading_data(salmon_data_dir)
