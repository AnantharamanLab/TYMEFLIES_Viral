#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    from collections import defaultdict
    from pathlib import Path
    import subprocess
    from subprocess import DEVNULL, STDOUT, check_call    
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: Find the most abundant AMG from viruses that infect Nanopelagicales


# Step 1 Get all the viruses that infect Nanopelagicales
Nanopelagicales_viruses = []
with open('Host_prediction/Viral_gn2host_tax_final.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        if 'Nanopelagicales' in tmp[1]:
            Nanopelagicales_viruses.append(tmp[0])    
            

# Step 2 Get all the virus representatives that infect Nanopelagicales
Species = {}  # Dictionary to store species information, gn_rep => gns
gn2species = {} # gn => species
with open('/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt', 'r') as file:
    for line in file:
        line = line.strip()
        tmp = line.split("\t")
        gn_rep = tmp[0]
        gns = tmp[1]
        Species[gn_rep] = gns
        gn_list = gns.split(',')
        for gn in gn_list:
            gn2species[gn] = gn_rep    

Nanopelagicales_virus_rep = []
for virus in Nanopelagicales_viruses:
    virus_rep = gn2species[virus]
    if virus_rep not in Nanopelagicales_virus_rep:
        Nanopelagicales_virus_rep.append(virus_rep)    


# Step 3 Get the most abundant AMG from these virus representatives
virus2amg_kos = defaultdict(list) # virus => [amg_kos]
with open('AMG_analysis/AMG_summary.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('Protein'):
            tmp = line.split('\t')
            virus = tmp[0].split('__Ga')[0]
            amg_ko = tmp[2]
            virus2amg_kos[virus].append(amg_ko)
            
amg_ko2num = defaultdict(int)
for virus in Nanopelagicales_virus_rep:
    if virus in virus2amg_kos:
        amg_kos = virus2amg_kos[virus]
        for amg_ko in amg_kos:
            amg_ko2num[amg_ko] += 1

print(amg_ko2num)           
    
    