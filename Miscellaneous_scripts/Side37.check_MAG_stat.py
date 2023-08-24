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
    
    
# Aim: Check the MAGs that passed the completeness and contamination criterion


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


# Step 2 Get the set of MAGs that passed the completeness and contamination criterion
MAG_set_1 = set() # MAGs that passed the completeness and contamination criterion
for MAG in MAG_stat:
    stat = MAG_stat[MAG]
    completeness, contamination = float(stat[7]), float(stat[8])
    if completeness >= 50 and contamination <= 10:
        MAG_set_1.add(MAG)
        

# Step 3 Get the set of MAGs that have num.in.cluster
MAG_set_2 = set() # MAGs that have num.in.cluster
for MAG in MAG_stat:
    stat = MAG_stat[MAG]
    num_in_cluster = stat[15]
    if num_in_cluster != 'NA':
        MAG_set_2.add(MAG) 
        
        
# Step 4 Compare MAG_set_1 and MAG_set_2
if MAG_set_1 == MAG_set_2:
    print("MAGs that passed the completeness and contamination criterion are the same with MAGs that have num.in.cluster")
    
print(f"The number of MAGs that passed the completeness and contamination criterion is {len(MAG_set_1)}")
print(f"The number of MAGs have num.in.cluster is {len(MAG_set_2)}")   
    
         
    


       