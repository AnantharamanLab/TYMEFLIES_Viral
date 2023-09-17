#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    from pathlib import Path 
    import subprocess
    from subprocess import DEVNULL, STDOUT, check_call 
    from collections import defaultdict  
    from datetime import datetime    
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)
    
# Aim: Parse to get ahbD-containing virus abundance from 0 day of Clearwater


# Step 1 Calculate the IMG2clearwater_day dict
## Step 1.1 Store each year's Clearwater state datetime
year2clearwater_start_day = {}
with open('season_start_dates.txt', 'r') as lines:
    for line in lines:
        if not line.startswith('Year'):
            line = line.rstrip('\n')
            year, clearwater_start_day = line.split('\t')[0], line.split('\t')[1]
            year2clearwater_start_day[year] = clearwater_start_day
           
## Step 1.2 Calculate to get the IMG2clearwater_day dict
### Get IMG2year, IMG2date, and year2IMGs dicts
IMG2date = {} # IMG => date
IMG2year = {} # IMG => year
year2IMGs = defaultdict(list) # year => [IMGs]
with open('/storage1/data11/TYMEFLIES_phage/TYMEFLIES_metagenome_info.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('IMG'):
            tmp = line.split('\t')
            IMG, date = tmp[0], tmp[8]
            IMG2date[IMG] = date
            year = date.split('-')[0]
            IMG2year[IMG] = year
            year2IMGs[year].append(IMG)
lines.close() 

IMG2clearwater_day = {} # IMG => clearwater_day
for IMG in IMG2date:
    date = IMG2date[IMG]
    year = IMG2year[IMG]
    clearwater_start_day = year2clearwater_start_day[year]
    
    date1 = datetime.strptime(clearwater_start_day, "%Y-%m-%d")
    date2 = datetime.strptime(date, "%Y-%m-%d")
    clearwater_day = date2 - date1
    IMG2clearwater_day[IMG] = clearwater_day.days
   
      
# Step 2 Get ahbD-containing viral gn list (species representatives)
Species = {}  # $gn_rep => $gns
AMG_summary = {}  # $pro => $ko
KOs = {}  # $ko => 1
AMG_containing_viral_gn = {}  # $gn => 1
Old_gene2new_gene_map = {}  # $gene_old => $gene_new
AhbD_containing_viral_gn = {}  # $gn => 1

## Step 2.1 Store species info
with open("/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt", "r") as file:
    for line in file:
        line = line.strip()
        tmp = line.split("\t")
        gn_rep = tmp[0]
        gns = tmp[1]
        Species[gn_rep] = gns

## Step 2.2 Store AMG KO information
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
            img_id = pro.split("_")[0]
            KOs[ko] = ko_detail
            gn = pro.split("__Ga")[0]
            AMG_containing_viral_gn[gn] = 1

### Change the old gene to new gene
with open("New_gene2old_gene_map.txt", "r") as file:
    for line in file:
        line = line.strip()
        tmp = line.split("\t")
        gene_new = tmp[0]
        gene_old = tmp[1]
        Old_gene2new_gene_map[gene_old] = gene_new

for pro in sorted(AMG_summary.keys()):
    if pro in Old_gene2new_gene_map:
        ko = AMG_summary[pro]
        gene_new = Old_gene2new_gene_map[pro]
        del AMG_summary[pro]
        AMG_summary[gene_new] = ko

## Step 2.3 Get ahbD-containing viral genome
for pro in sorted(AMG_summary.keys()):
    gn = pro.split("__Ga")[0]
    ko = AMG_summary[pro]
    if gn in Species and ko == "K22227":
        AhbD_containing_viral_gn[gn] = 1
        
       
# Step 3 Store viral_gn2IMG2cov dict
viral_gn2IMG2cov = defaultdict(dict)
Header = []  # Store the header line
with open('MetaPop/Viral_gn2IMG2cov_norm_filtered.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Head'):
            Header = line.split('\t')
        else:
            tmp = line.split('\t')
            for i in range(1, len(tmp)):
                viral_gn = tmp[0]
                IMG = Header[i]
                cov = tmp[i]
                if cov == 'NA':
                    cov = '0'
                viral_gn2IMG2cov[viral_gn][IMG] = float(cov)


# Step 4 Calculate IMG2ahbD_containing_viral_gn_cov dict
IMG2ahbD_containing_viral_gn_cov = {} # IMG => ahbD_containing_viral_gn_cov (The sum cov of ahbD-containing viral gn (rep))            
for IMG in IMG2date:
    ahbD_containing_viral_gn_cov = 0
    for viral_gn in viral_gn2IMG2cov:
        if viral_gn in AhbD_containing_viral_gn:
            cov = viral_gn2IMG2cov[viral_gn][IMG]
            ahbD_containing_viral_gn_cov += float(cov)
    IMG2ahbD_containing_viral_gn_cov[IMG] = ahbD_containing_viral_gn_cov    
         
        
# Step 5 Calculate clearwater_day2ahbD_containing_viral_gn_cov dicts for each year
os.mkdir("virus_n_MAG_tax_association/ahbD_containing_virus")
for year in year2IMGs:
    IMGs = year2IMGs[year]
    clearwater_day2ahbD_containing_viral_gn_cov = {} # int(clearwater_day) => ahbD_containing_viral_gn_cov
    for IMG in IMGs:
        clearwater_day = int(IMG2clearwater_day[IMG])
        ahbD_containing_viral_gn_cov = IMG2ahbD_containing_viral_gn_cov[IMG]
        clearwater_day2ahbD_containing_viral_gn_cov[clearwater_day] = ahbD_containing_viral_gn_cov

    ## Write down the result
    f = open(f"virus_n_MAG_tax_association/ahbD_containing_virus/{year}.ahbD_containing_viral_gn_cov.txt", 'w')
    f.write('clearwater_day\tahbD_containing_viral_gn_cov\n')
    for clearwater_day in sorted(clearwater_day2ahbD_containing_viral_gn_cov.keys()):
        line = str(clearwater_day) + '\t' + str(clearwater_day2ahbD_containing_viral_gn_cov[clearwater_day])
        f.write(line + '\n')