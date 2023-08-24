#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    import statistics
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: Get the binned scaffold ratio ranges from 465 metagenomic datasets


# Step 1 Get all the folders contain viral genomes
all_viral_gn_folder_addrs = glob('/storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed')


# Step 2 Get the binned scaffold ratios from each folder
binned_scaffold_ratios = []
for each_viral_gn_folder_addr in all_viral_gn_folder_addrs:
    binned_scaffold_no = 0
    unbinned_scaffold_no = 0
    all_fasta_addrs = glob(f'{each_viral_gn_folder_addr}/*.fasta')
    for each_fasta_addr in all_fasta_addrs:
        if 'vRhyme_unbinned' in each_fasta_addr: # for unbinned scaffolds
            unbinned_scaffold_no = unbinned_scaffold_no + 1
        else: # for binned scaffolds
            with open(each_fasta_addr, 'r') as lines:
                for line in lines:
                    line = line.rstrip('\n')
                    if line.startswith('>'):
                        binned_scaffold_no = binned_scaffold_no + 1
            lines.close()
    binned_scaffold_ratio = binned_scaffold_no / (binned_scaffold_no + unbinned_scaffold_no)
    binned_scaffold_ratio = round(binned_scaffold_ratio, 3) # Get only three digitals
    binned_scaffold_ratios.append(binned_scaffold_ratio)
    
    
# Step 3 Get the range, mean, median of the list of binned_scaffold_ratios
## Calculate the range as percentages
range_str = f"{round(min(binned_scaffold_ratios)*100, 1)}% - {round(max(binned_scaffold_ratios)*100, 1)}%"

## Calculate the mean as a percentage
mean = round(sum(binned_scaffold_ratios) / len(binned_scaffold_ratios) * 100, 1)

## Calculate the median as a percentage
median = round(statistics.median(binned_scaffold_ratios) * 100, 1)

## Print the result
print("Binned scaffold ratios range:", range_str)
print("Binned scaffold ratios mean:", mean, "%")
print("Binned scaffold ratios median:", median, "%")                    