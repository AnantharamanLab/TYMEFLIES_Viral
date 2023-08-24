#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    import re
    from glob import glob
    warnings.filterwarnings('ignore')
    from collections import defaultdict
    from pathlib import Path
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: (1) Get the stat of MAGs that passed the CheckM completess and contamination criterion
#      (2) Get the stat of TYMEFLIES rep MAGs

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
    

# Step 1 Store the information from Robin_MAG_stat.txt, incude tax and scaffolds
TYMEFLIES_MAG_stat = {}  # mag => [lineage, completeness, contamination, scaffolds]
TYMEFLIES_rep_MAG = [] # Store the list of TYMEFLIES_rep_MAG
with open('/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.txt', 'r') as file:
    for line in file:
        line = line.strip()
        if not line.startswith('tymeflies'):
            tmp = line.split('\t')
            mag = tmp[5]
            winner = tmp[14]
            num_in_cluster = tmp[15]            
            completeness, contamination = tmp[7], tmp[8]
            if num_in_cluster != "NA":
                img = re.search(r'_(33\d+?)_', mag).group(1)
                Contigs = []  # Store all the contigs into an array
                MAG_addr = "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/" + img + "/" + mag + ".fasta"
                MAG_seq = store_seq(MAG_addr)  # Assuming store_seq function exists
                for header in sorted(MAG_seq.keys()):
                    contig = header[1:]
                    Contigs.append(contig)
                lineage = ';'.join(tmp[16:23])
                scaffolds = ','.join(Contigs)  # Store all scaffolds in each MAG
                TYMEFLIES_MAG_stat[mag] = [lineage, completeness, contamination, scaffolds]
                if winner == 'TRUE':
                    TYMEFLIES_rep_MAG.append(mag)
      
'''     
# Step 2 Write down the Robin_MAG_stat.CheckM_passed.txt           
output_file = 'Robin_MAG_stat.CheckM_passed.txt'
header_line = "IMG Bin ID\tGTDB-TK lineage\tBin Completeness\tBin Contamination\tScaffolds\n"

with open(output_file, 'w') as file:
    file.write(header_line)
    for mag, info in TYMEFLIES_MAG_stat.items():
        lineage, completeness, contamination, scaffolds = info
        line = f"{mag}\t{lineage}\t{completeness}\t{contamination}\t{scaffolds}\n"
        file.write(line)   
'''

# Step 3 Write down the Robin_MAG_stat.rep_MAG.txt          
output_file = 'Robin_MAG_stat.rep_MAG.txt'
header_line = "IMG Bin ID\tGTDB-TK lineage\tBin Completeness\tBin Contamination\tScaffolds\n"

with open(output_file, 'w') as file:
    file.write(header_line)
    for mag, info in TYMEFLIES_MAG_stat.items():
        if mag in TYMEFLIES_rep_MAG:
            lineage, completeness, contamination, scaffolds = info
            line = f"{mag}\t{lineage}\t{completeness}\t{contamination}\t{scaffolds}\n"
            file.write(line)         