import csv
from pathlib import Path
from pympler.asizeof import asizeof


def salmon_reading_data(data_dir: str) -> dict:
	# Finding data files
	data_paths = Path(data_dir).rglob('quant.genes.sf')
	data_dict = {}
	
	# Reading data from files and combining them into one big object
	for path in data_paths:
		sample_name = path.parts[-2]
		reader = csv.reader(open(path), delimiter='\t')
		first_line = True
		name_index = 0
		counts_index = 0
		data = []
		for line in reader:
			if first_line:
				name_index = line.index('Name')
				counts_index = line.index('NumReads')
				first_line = False
				continue
			gene_name = line[name_index]
			gene_count = float(line[counts_index])
			data.append((gene_name, gene_count))
		data_dict[sample_name] = data
	return data_dict


# Directory of data, later to be changed for agrpass
salmon_data_dir = 'SALMON_OUT'

data_dict = salmon_reading_data(salmon_data_dir)

# Prints size for testing purposes
size = asizeof(data_dict)
print(size/1000**2, 'Mb')
