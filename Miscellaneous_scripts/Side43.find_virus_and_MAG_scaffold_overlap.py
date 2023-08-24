#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    from collections import defaultdict
    import subprocess
    from subprocess import DEVNULL, STDOUT, check_call    
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: Find the virus and MAG scaffold overlap rate


# Step 1 Get all the viral rep gn scaffold list
viral_rep_gn_scaffold_list = []
with open('/storage1/data11/TYMEFLIES_phage/reference_fasta_for_metapop/All_phage_species_rep_gn_containing_AMG.fasta', 'r') as lines:
    for line in lines:
        if line.startswith('>'):
            line = line.rstrip('\n')
            scaffold = line.split('__')[2]
            if 'fragment' in scaffold:
                scaffold = scaffold.split('_frag')[0]
            viral_rep_gn_scaffold_list.append(scaffold)    

# Step 2 Get all the rep MAG scaffold list   
rep_MAG_scaffold_list = []
with open('TYMEFLIES_rep_MAG.fasta', 'r') as lines:
    for line in lines:
        if line.startswith('>'):
            line = line.rstrip('\n').replace('>', '',1)
            rep_MAG_scaffold_list.append(line)

# Step 3 Find the overlaped items and report numbers

# Find the overlapping item numbers
overlapping_items = set(viral_rep_gn_scaffold_list) & set(rep_MAG_scaffold_list)

# Report the item numbers of all three lists
print("viral_rep_gn_scaffold_list:", len(viral_rep_gn_scaffold_list))
print("rep_MAG_scaffold_list:", len(rep_MAG_scaffold_list))
print("Overlapping items:", len(overlapping_items))
    
    