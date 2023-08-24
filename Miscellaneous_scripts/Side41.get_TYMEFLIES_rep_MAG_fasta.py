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
    
    
# Aim: Get the concatenated fasta file for all TYMEFLIES representative MAGs


# Step 1 Get the rep_MAG_list
rep_MAG_list = []
with open('/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.txt','r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('tymeflies'):
            tmp = line.split('\t')
            MAG, winner = tmp[5], tmp[14]
            if winner == 'TRUE':
                rep_MAG_list.append(MAG)

print(f"The total number of rep MAGs is {len(rep_MAG_list)}")          

# Step 2 Concatenate all fasta files
MAG_addrs = []
for MAG in rep_MAG_list:
    MAG_addr = f"/storage1/data11/TYMEFLIES_phage/Robin_MAGs/All_passed_MAGs/{MAG}.fasta"
    if os.path.exists(MAG_addr):
        MAG_addrs.append(MAG_addr)        

output_file = f"/storage1/data11/TYMEFLIES_phage/TYMEFLIES_rep_MAG.fasta"

with open(output_file, 'w') as outfile:
    for fasta_file in MAG_addrs:
        with open(fasta_file, 'r') as infile:
            for line in infile:
                outfile.write(line)
