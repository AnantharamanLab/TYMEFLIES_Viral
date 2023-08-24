#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    from collections import defaultdict
    from pathlib import Path
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: Get AMG proteins from six IMG results


def store_seq(input_seq_file): # The input sequence file should be a file with full path
    head = "" # Store the header line
    seq_dict = {} # Store the sequence dict
    
    with open(input_seq_file, "r") as seq_lines:
        for line in seq_lines:
            line = line.rstrip("\n") # Remove "\n" in the end
            if ">" in line:
                head = line.split(None, 1)[0] # Cut at the first " " or "\t", use the first part
                seq_dict[head] = ""                 
            else:
                seq_dict[head] += line
            
    seq_lines.close()
    
    return seq_dict
    

# Step 1 Get all AMG protein list from six IMG results
## Step 1.1 Store the six IMG ID
IMGID = {}
IMGID["3300049597"] = 1
IMGID["3300049596"] = 1
IMGID["3300048593"] = 1
IMGID["3300049595"] = 1
IMGID["3300049594"] = 1
IMGID["3300049629"] = 1

## Step 1.2 Open the file containing the AMG list and read the AMGs into a list
All_AMG_list = [] # Store all the AMG list
with open('AMG_analysis/AMG_summary.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('protein'):
            AMG = line.split('\t')[0]
            AMG_ID = AMG.split('__', 1)[0]
            if AMG_ID in IMGID:
                All_AMG_list.append(AMG)
lines.close()


# Step 2 Get AMG proteins
## Step 2.1 Store all viral proteins
all_virus_protein_seq = {} # The seq dict for all virus proteins
all_viral_faa_addrs = glob('/storage1/data11/TYMEFLIES_phage/*/vRhyme_best_bins_fasta_parsed/*.faa')
for each_viral_faa_addr in all_viral_faa_addrs:
    if each_viral_faa_addr.split('/')[4] in IMGID:
        each_viral_faa_seq = store_seq(each_viral_faa_addr)
        all_virus_protein_seq.update(each_viral_faa_seq)
        
## Step 2.2 Get AMG protein
all_virus_AMG_protein_seq = {} # The seq dict for all virus AMG proteins
for header in all_virus_protein_seq:
    header_wo_array = header.replace('>', '', 1)
    if header_wo_array in All_AMG_list:
        all_virus_AMG_protein_seq[header] = all_virus_protein_seq[header]
        
## Step 2.3 Write down the AMG proteins
f = open('AMG_protein_from_6_IMG.faa', 'w')
for header in all_virus_AMG_protein_seq:
    f.write(f'{header}\n{all_virus_AMG_protein_seq[header]}\n')
f.close()    

                







       