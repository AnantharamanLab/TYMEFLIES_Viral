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
    
# Aim: Parse to get psbA virus and MAG abundance from 0 day of Early Summer


# Step 1 Calculate the IMG2earlysummer_day dict
## Step 1.1 Store each year's Early Summer start datetime
year2earlysummer_start_day = {}
with open('season_start_dates.txt', 'r') as lines:
    for line in lines:
        if not line.startswith('Year'):
            line = line.rstrip('\n')
            year, earlysummer_start_day = line.split('\t')[0], line.split('\t')[3]
            year2earlysummer_start_day[year] = earlysummer_start_day
           
## Step 1.2 Calculate to get the IMG2earlysummer_day dict
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

IMG2earlysummer_day = {} # IMG => earlysummer_day
for IMG in IMG2date:
    date = IMG2date[IMG]
    year = IMG2year[IMG]
    earlysummer_start_day = year2earlysummer_start_day[year]
    
    date1 = datetime.strptime(earlysummer_start_day, "%Y-%m-%d")
    date2 = datetime.strptime(date, "%Y-%m-%d")
    earlysummer_day = date2 - date1
    IMG2earlysummer_day[IMG] = earlysummer_day.days
    
    
# Step 2 Calculate IMG2Cyanobiaceae_MAG_cov dict
## Step 2.1 Store MAG2IMG2abun dict
MAG2IMG2abun = defaultdict(dict)
Header2 = []  # Store the header line
with open('MAG_abundance/MAG2IMG2abun.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Head'):
            Header2 = line.split('\t')
        else:
            tmp = line.split('\t')
            for i in range(1, len(tmp)):
                MAG = tmp[0]
                IMG = Header2[i]
                abun = tmp[i]
                MAG2IMG2abun[MAG][IMG] = abun
                
## Step 2.2 Store Cyanobiaceae_MAG list
Cyanobiaceae_MAG = []
with open('/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.CheckM_passed.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('IMG'):
            tmp = line.split('\t')
            MAG, tax = tmp[0], tmp[1]
            if 'f__Cyanobiaceae' in tax:
                Cyanobiaceae_MAG.append(MAG)
        
## Step 2.3  Calculate IMG2Cyanobiaceae_MAG_cov dict          
IMG2Cyanobiaceae_MAG_cov = {} # IMG => Cyanobiaceae_MAG_cov (The sum cov of Cyanobiaceae MAGs)            
for IMG in IMG2date:
    Cyanobiaceae_MAG_cov = 0
    for MAG in Cyanobiaceae_MAG:
        cov = MAG2IMG2abun[MAG][IMG]
        Cyanobiaceae_MAG_cov += float(cov)
    IMG2Cyanobiaceae_MAG_cov[IMG] = Cyanobiaceae_MAG_cov


# Step 3 Calculate earlysummer_day2Cyanobiaceae_MAG_cov dict for each year 
#os.mkdir('virus_n_MAG_tax_association')
os.mkdir('virus_n_MAG_tax_association/psbA_virus_n_MAG')
for year in year2IMGs:
    IMGs = year2IMGs[year]
    earlysummer_day2Cyanobiaceae_MAG_cov = {} # int(earlysummer_day) => Cyanobiaceae_MAG_cov
    for IMG in IMGs:
        earlysummer_day = int(IMG2earlysummer_day[IMG])
        Cyanobiaceae_MAG_cov = IMG2Cyanobiaceae_MAG_cov[IMG]
        earlysummer_day2Cyanobiaceae_MAG_cov[earlysummer_day] = Cyanobiaceae_MAG_cov

    ## Write down the result
    f = open(f"virus_n_MAG_tax_association/psbA_virus_n_MAG/{year}.Cyanobiaceae_MAG_cov.txt", 'w')
    f.write('earlysummer_day\tCyanobiaceae_MAG_cov\n')
    for earlysummer_day in sorted(earlysummer_day2Cyanobiaceae_MAG_cov.keys()):
        line = str(earlysummer_day) + '\t' + str(earlysummer_day2Cyanobiaceae_MAG_cov[earlysummer_day])
        f.write(line + '\n') 
    f.close()    
        
# Step 4 Get psbA-containing viral gn list (species representatives)
Species = {}  # $gn_rep => $gns
AMG_summary = {}  # $pro => $ko
KOs = {}  # $ko => 1
AMG_containing_viral_gn = {}  # $gn => 1
Old_gene2new_gene_map = {}  # $gene_old => $gene_new
PsbA_containing_viral_gn = {}  # $gn => 1

## Step 4.1 Store species info
with open("/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt", "r") as file:
    for line in file:
        line = line.strip()
        tmp = line.split("\t")
        gn_rep = tmp[0]
        gns = tmp[1]
        Species[gn_rep] = gns

## Step 4.2 Store AMG KO information
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

## Step 4.3 Get psbA-containing viral genome
for pro in sorted(AMG_summary.keys()):
    gn = pro.split("__Ga")[0]
    ko = AMG_summary[pro]
    if gn in Species and ko == "K02703":
        PsbA_containing_viral_gn[gn] = 1
        
       
# Step 5 Store viral_gn2IMG2cov dict
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
                viral_gn2IMG2cov[viral_gn][IMG] = cov


# Step 6 Calculate IMG2psbA_viral_gn_cov dict
IMG2psbA_viral_gn_cov = {} # IMG => psbA_viral_gn_cov (The sum cov of psbA-containing viral gn (rep))            
for IMG in IMG2date:
    psbA_viral_gn_cov = 0
    for viral_gn in viral_gn2IMG2cov:
        if viral_gn in PsbA_containing_viral_gn:
            cov = viral_gn2IMG2cov[viral_gn][IMG]
            psbA_viral_gn_cov += float(cov)
    IMG2psbA_viral_gn_cov[IMG] = psbA_viral_gn_cov    
         
        
# Step 7 Calculate earlysummer_day2psbA_viral_gn_cov dicts for each year
for year in year2IMGs:
    IMGs = year2IMGs[year]
    earlysummer_day2psbA_viral_gn_cov = {} # int(earlysummer_day) => psbA_viral_gn_cov
    for IMG in IMGs:
        earlysummer_day = int(IMG2earlysummer_day[IMG])
        psbA_viral_gn_cov = IMG2psbA_viral_gn_cov[IMG]
        earlysummer_day2psbA_viral_gn_cov[earlysummer_day] = psbA_viral_gn_cov

    ## Write down the result
    f = open(f"virus_n_MAG_tax_association/psbA_virus_n_MAG/{year}.psbA_viral_gn_cov.txt", 'w')
    f.write('earlysummer_day\tpsbA_viral_gn_cov\n')
    for earlysummer_day in sorted(earlysummer_day2psbA_viral_gn_cov.keys()):
        line = str(earlysummer_day) + '\t' + str(earlysummer_day2psbA_viral_gn_cov[earlysummer_day])
        f.write(line + '\n')
    f.close()    
        
        
# Step 8 Store non-psbA containing virus to IMG 2 cov norm filtered
no_AMG_viral_gn2IMG2cov_norm_filtered = defaultdict(dict) # viral_gn => IMG => cov_norm_filtered
Header3 = [] # Store the header line
with open('MetaPop/no_AMG_viral_gn2IMG2cov_norm_filtered.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Head'):
            tmp = line.split('\t')
            Header3 = tmp
        else:
            tmp = line.split('\t')
            viral_gn = tmp[0]
            for i in range(1, len(tmp)):
                IMG = Header3[i]
                cov_norm_filtered = tmp[i]
                no_AMG_viral_gn2IMG2cov_norm_filtered[viral_gn][IMG] = float(cov_norm_filtered)   


# Step 9 Store species to host family (taxonomy) dict
# Step 9.1 Store viral_gn2host_family dict and viral_gn2Cyanobiaceae dict
viral_gn2host_family = {} # viral_gn => host_family
viral_gn2Cyanobiaceae = set() # Set to store viral_gn with host as Cyanobiaceae
with open('/storage1/data11/TYMEFLIES_phage/Host_prediction/Viral_gn2host_tax_final.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        viral_gn, host_family = tmp[0], tmp[1]
        host_family = host_family.split(';g__', 1)[0]
        viral_gn2host_family[viral_gn] = host_family
        
        # Check if the host is Cyanobiaceae and add the viral_gn to the set
        if host_family == 'd__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__PCC-6307;f__Cyanobiaceae':
            viral_gn2Cyanobiaceae.add(viral_gn)    
            
        
## Step 9.2 Store non-psbA containing species with Cyanobiaceae host
non_psbA_species_with_Cyanobiaceae_host = set()
for viral_gn in viral_gn2Cyanobiaceae:
    if viral_gn in no_AMG_viral_gn2IMG2cov_norm_filtered:
        non_psbA_species_with_Cyanobiaceae_host.add(viral_gn)
    elif viral_gn in viral_gn2IMG2cov and viral_gn not in PsbA_containing_viral_gn:
        non_psbA_species_with_Cyanobiaceae_host.add(viral_gn)


# Step 10 Calculate IMG2non_psbA_Cyanobiaceae_virus_cov dict
IMG2non_psbA_Cyanobiaceae_virus_cov = {} # IMG => non_psbA_Cyanobiaceae_virus_cov          
for IMG in IMG2date:
    non_psbA_Cyanobiaceae_virus_cov = 0
    for viral_gn in no_AMG_viral_gn2IMG2cov_norm_filtered:
        if viral_gn in non_psbA_species_with_Cyanobiaceae_host:
            cov = no_AMG_viral_gn2IMG2cov_norm_filtered[viral_gn][IMG]
            non_psbA_Cyanobiaceae_virus_cov += cov
    for viral_gn in viral_gn2IMG2cov:
        if viral_gn in non_psbA_species_with_Cyanobiaceae_host:
            cov = viral_gn2IMG2cov[viral_gn][IMG]
            non_psbA_Cyanobiaceae_virus_cov += cov            
    IMG2non_psbA_Cyanobiaceae_virus_cov[IMG] = non_psbA_Cyanobiaceae_virus_cov
    
  
# Step 11 Calculate earlysummer_day2non_psbA_Cyanobiaceae_virus_cov dicts for each year
for year in year2IMGs:
    IMGs = year2IMGs[year]
    earlysummer_day2non_psbA_Cyanobiaceae_virus_cov = {} # int(earlysummer_day) => non_psbA_Cyanobiaceae_virus_cov
    for IMG in IMGs:
        earlysummer_day = int(IMG2earlysummer_day[IMG])
        non_psbA_Cyanobiaceae_virus_cov = IMG2non_psbA_Cyanobiaceae_virus_cov[IMG]
        earlysummer_day2non_psbA_Cyanobiaceae_virus_cov[earlysummer_day] = non_psbA_Cyanobiaceae_virus_cov

    ## Write down the result
    f = open(f"virus_n_MAG_tax_association/psbA_virus_n_MAG/{year}.non_psbA_Cyanobiaceae_virus_cov.txt", 'w')
    f.write('earlysummer_day\tnon_psbA_Cyanobiaceae_virus_cov\n')
    for earlysummer_day in sorted(earlysummer_day2non_psbA_Cyanobiaceae_virus_cov.keys()):
        line = str(earlysummer_day) + '\t' + str(earlysummer_day2non_psbA_Cyanobiaceae_virus_cov[earlysummer_day])
        f.write(line + '\n')  
    f.close()    
