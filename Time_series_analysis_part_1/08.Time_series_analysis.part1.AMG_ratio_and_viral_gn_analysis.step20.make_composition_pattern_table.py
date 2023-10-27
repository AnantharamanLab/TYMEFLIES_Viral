#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    import datetime
    from collections import defaultdict
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: Make the composition pattern table based on viral gn abundance at the family level


# Step 1 Get the viral_gn2IMG2cov dict
viral_gn2IMG2cov = defaultdict(dict)
IMG_set = set()
Header = [] # Store the header line
with open('MetaPop/Viral_gn2IMG2cov_norm_filtered.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Head'):
            Header = line.split('\t')    
        else:
            tmp = line.split('\t')
            for i in range(1, len(tmp)):
                IMG = Header[i]
                viral_gn = tmp[0]
                cov = tmp[i]
                if cov == 'NA':
                    cov = 0
                viral_gn2IMG2cov[viral_gn][IMG] = cov 
                IMG_set.add(IMG)
                
Header2 = [] # Store the header line
with open('MetaPop/no_AMG_viral_gn2IMG2cov_norm_filtered.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Head'):
            Header = line.split('\t')    
        else:
            tmp = line.split('\t')
            for i in range(1, len(tmp)):
                IMG = Header[i]
                viral_gn = tmp[0]
                cov = tmp[i]
                if cov == 'NA':
                    cov = 0
                viral_gn2IMG2cov[viral_gn][IMG] = cov                 


# Step 2 Get family2IMG2abun dict
## Step 2.1 Store family2viral_gns dict
family2viral_gns = {} # family => [viral_gns]
j = 1  # The number of cluster
with open('Cluster_phage_genomes/family_clusters.txt', 'r') as file:
    for line in file:
        line = line.rstrip('\n')
        family = f"Family_{j:05d}"
        family2viral_gns[family] = line.split('\t')
        j += 1
        
## Step 2.2 Get family2IMG2abun dict
family2IMG2abun = defaultdict(dict)
for family in family2viral_gns:
    for IMG in IMG_set:
        viral_gns = family2viral_gns[family]
        for viral_gn in viral_gns:
            cov = 0
            if viral_gn2IMG2cov[viral_gn].get(IMG):
                cov = float(viral_gn2IMG2cov[viral_gn][IMG])
            family2IMG2abun[family][IMG] = family2IMG2abun[family].get(IMG, 0) + cov
        
 
# Step 3 Generate the compositional pattern table
## Step 3.1 Get IMG2season dict
IMG2season = {} # IMG => season
with open('TYMEFLIES_metagenome_info.txt', 'r') as lines:
    for line in lines:
        if not line.startswith('IMG'):
            line = line.rstrip('\n')
            tmp = line.split('\t')
            IMG2season[tmp[0]] = tmp[10]
            
## Step 3.2 Write down the compositional_pattern_table.txt
with open('compositional_pattern_table.txt', 'w') as file:
    ## Write table header
    header = '\t'.join(sorted(family2viral_gns.keys()))
    file.write('head\tseason\t' + header + '\n')

    ## Write table rows
    for IMG in sorted(IMG2season.keys()):
        row_values = [str(family2IMG2abun[family].get(IMG, '')) for family in sorted(family2IMG2abun.keys())]
        row = IMG + '\t' + IMG2season[IMG] + '\t' + '\t'.join(row_values) + '\n'
        file.write(row)             
