# Diff-exp-combine-harvester
This program is for robust differential expression analysis with single terminal command

# Description
The program uses R packages that calculate differential expression from output by either Salmon or RSEM.
After that initial calculations will be done, and it will produce some graphs from inside R, it will compare what genes were output by both packages.
It will merge information from both scripts into one tsv file output where user wants it, but also will ask Biomart service about what function was assigned to each gene.
The functions and their will be plotted and saved to tsv file.

# Quick start
To use this program, the best choice is to clone repo to chosen directory and use -h to get some help and all parameters.\
\
There is also an option to use example data provided with this program. There are two datasets (one from RSEM and one from Salmon) and you can choose either one
```
python main.py Example_data/SALMON_OUT salmon Example_data/salmon_experiment_design.txt
```

## Options
All options available:
- dirname - path to directory containing results of one of supported programs
- input_type - one of (salmon, rsem). Specifies what program was used to calculate gene counts
- design_filename - path to file explaining experiment design, see EXPERIMENT DESIGN FILE section on how to write one
- -o/--out_dir - Optional argument used to specify what directory save results to, defaults to cwd, will be created if non-existent
- -fc/--fc_cut_off - Cut off for log fc values when selecting genes to be analysed by comparison, defaults to 1.5

## Experiment design file
This is file format I came up with to easily tell program what is the structure of experiment. It can be easily edited and can be made by other scripts if needed.
Can also be used to limit what samples will be used. It can only accept two groups at once, but there is no limit to number of samples or names of groups, those will be used to name output tsv files. \
\
Example file
``` 
reference: SRR1747399, SRR1747397, SRR1747395
treatment: SRR1747301, SRR1747299, SRR1747303
```

### Notes
This is a repo of my project done during my apprenticeship in [data2biology](https://data2biology.com/)
