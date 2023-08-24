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
    
    
# Aim: Check the non-psbA virus and psbA virus (species rep) that can infect "f__Cyanobiaceae"


# Step 1 Store species info
Species = {} # gn_rep => gns (for example, "3300033816__vRhyme_113,3300034101__vRhyme_unbinned4136,3300042863__vRhyme_unbinned973")
with open("/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt", "r") as file:
    for line in file:
        line = line.strip()
        tmp = line.split("\t")
        gn_rep = tmp[0]
        gns = tmp[1]
        Species[gn_rep] = gns

# Step 2 Get the non-psbA species rep and psbA species that can infect "f__Cyanobiaceae" numbers
## Step 2.1 Store AMG KO information
AMG_summary = {}
KOs = {}
IMG2date = {}
AMG_containing_viral_gn = defaultdict(list) # gn => [kos]
with open("AMG_analysis/AMG_summary.txt", "r") as file:
    for line in file:
        line = line.strip()
        if not line.startswith("Pro"):
            tmp = line.split("\t")
            pro = tmp[0]
            ko = tmp[2]
            ko_detail = tmp[3]
            date_n_season = tmp[1]
            AMG_summary[pro] = ko
            img_id = pro.split("__")[0]
            IMG2date[img_id] = date_n_season
            KOs[ko] = ko_detail

            gn = pro.split("__Ga")[0]
            AMG_containing_viral_gn[gn].append(ko)
            
## Step 2.2 Store viral_gn2host_tax dict
viral_gn2host_tax = {} # viral_gn => host_tax   
with open('Host_prediction/Viral_gn2host_tax_final.txt', 'r') as lines:
     for line in lines:
         line = line.rstrip('\n')
         tmp = line.split('\t')
         viral_gn2host_tax[tmp[0]] = tmp[1]

## Step 2.3 Store non_psbA_species_rep and psbA_species_rep that can infect "f__Cyanobiaceae"
non_psbA_species_rep = []
psbA_species_rep = []
for species_rep in Species:
    if species_rep in viral_gn2host_tax and 'f__Cyanobiaceae' in viral_gn2host_tax[species_rep]:
        # Store psbA-containing species rep that can infect "f__Cyanobiaceae"
        if species_rep in AMG_containing_viral_gn and 'K02703' in AMG_containing_viral_gn[species_rep]:
            psbA_species_rep.append(species_rep)  
        # Store non psbA-containing species rep that can infect "f__Cyanobiaceae"
        else:
            non_psbA_species_rep.append(species_rep)  
            
print(f"The number of psbA-containing species rep that can infect Cyanobiaceae is : {len(psbA_species_rep)}")            
print(f"The number of non psbA-containing species rep that can infect Cyanobiaceae is : {len(non_psbA_species_rep)}")      
    