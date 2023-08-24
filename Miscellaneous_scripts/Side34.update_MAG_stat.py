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
    
    
# Aim: Update the MAG stat table ('Robin_MAG_stat.txt') with the result of GTDB-Tk


# Step 1 Store the MAG stat table 
MAG_stat = {} # MAG => [whole line]; for example, ME2013-02-02s1D0_3300042325_group4_bin39 => [whole line]
header = '' # Store the header line of 'Robin_MAG_stat.txt'
with open('Robin_MAG_stat.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('tymeflies'):
            header = line
        else:
            tmp = line.split('\t')
            MAG = tmp[5]
            MAG_stat[MAG] = tmp
lines.close()   
         
    
# Step 2 Store the GTDB-Tk result
MAG2tax = {} # MAG => tax
## Step 2.1 Store the taxonomical result of archaeal MAGs
with open('all_MAGs_fasta_gtdbtk_output_dir/gtdbtk.ar53.summary.tsv', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('user_genome'):
            tmp = line.split('\t')
            MAG2tax[tmp[0]] = tmp[1]
lines.close()

## Step 2.2 Store the taxonomical result of bacterial MAGs
with open('all_MAGs_fasta_gtdbtk_output_dir/gtdbtk.bac120.summary.tsv', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('user_genome'):
            tmp = line.split('\t')
            MAG2tax[tmp[0]] = tmp[1]
lines.close()

           
# Step 3 Make new MAG stat table 
## Step 3.1 Replace the MAG_stat dict
for MAG in MAG_stat:    
    tax_list = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    if MAG in MAG2tax:
        tax = MAG2tax[MAG]
        if tax == 'Unclassified Archaea':
            tax_list = ['d__Archaea', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
        elif tax == 'Unclassified Bacteria':
            tax_list = ['d__Bacteria', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
        elif tax == 'Unclassified':
            tax_list = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
        else:
            tax_list = tax.split(';')
    
    stat = MAG_stat[MAG]     
    stat[16], stat[17], stat[18], stat[19], stat[20], stat[21], stat[22] = tax_list[0], tax_list[1], tax_list[2], tax_list[3], tax_list[4], tax_list[5], tax_list[6]
    MAG_stat[MAG] = stat
    
## Step 3.2 Write down the new MAG_stat file
f = open('Robin_MAG_stat.new.txt', 'w')
f.write(header + '\n')
for MAG in MAG_stat:
    line = '\t'.join(MAG_stat[MAG]) + '\n'
    f.write(line) 
f.close()    

       