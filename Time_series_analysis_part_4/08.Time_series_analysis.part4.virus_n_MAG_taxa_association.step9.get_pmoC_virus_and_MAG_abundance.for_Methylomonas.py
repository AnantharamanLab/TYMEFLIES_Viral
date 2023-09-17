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
    
# Aim: Parse to get Methylomonas virus (both pmoC-containing and no-pmoC-containing) 
#      and Methylomonas MAG abundance from 0 day of Late Summer


# Step 1 Calculate the IMG2latesummer_day dict
## Step 1.1 Store each year's Late Summer state datetime
year2latesummer_start_day = {}
with open('season_start_dates.txt', 'r') as lines:
    for line in lines:
        if not line.startswith('Year'):
            line = line.rstrip('\n')
            year, latesummer_start_day = line.split('\t')[0], line.split('\t')[4]
            year2latesummer_start_day[year] = latesummer_start_day
           
## Step 1.2 Calculate to get the IMG2latesummer_day dict
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

IMG2latesummer_day = {} # IMG => latesummer_day
for IMG in IMG2date:
    date = IMG2date[IMG]
    year = IMG2year[IMG]
    latesummer_start_day = year2latesummer_start_day[year]
    
    date1 = datetime.strptime(latesummer_start_day, "%Y-%m-%d")
    date2 = datetime.strptime(date, "%Y-%m-%d")
    latesummer_day = date2 - date1
    IMG2latesummer_day[IMG] = latesummer_day.days
    
    
# Step 2 Calculate Methylomonas MAG abundance
## Step 2.1 Store MAG2IMG2abun dict
MAG2IMG2abun = defaultdict(dict)
Header = []  # Store the header line
with open('virus_n_MAG_tax_association/MAG2IMG2abun.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('head'):
            Header = line.split('\t')
        else:
            tmp = line.split('\t')
            for i in range(1, len(tmp)):
                MAG = tmp[0]
                IMG = Header[i]
                abun = tmp[i]
                MAG2IMG2abun[MAG][IMG] = float(abun)
                
## Step 2.2 Store Methylomonas_MAG list
Methylomonas_genus_list = [ 
    'd__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Methylococcales;f__Methylomonadaceae;g__Methylomonas'
]

Methylomonas_MAG = []
with open('/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.rep_MAG.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('IMG'):
            tmp = line.split('\t')
            MAG, tax = tmp[0], tmp[1]
            genus = tax.split(';s__', 1)[0]
            if genus in Methylomonas_genus_list:
                Methylomonas_MAG.append(MAG)
        
## Step 2.3  Calculate IMG2Methylomonas_MAG_cov dict         
IMG2Methylomonas_MAG_cov = {} # IMG => Methylomonas_MAG_cov (The sum cov of Methylomonas MAGs)            
for IMG in IMG2date:
    Methylomonas_MAG_cov = 0
    for MAG in Methylomonas_MAG:
        cov = MAG2IMG2abun[MAG][IMG]
        Methylomonas_MAG_cov += cov
    IMG2Methylomonas_MAG_cov[IMG] = Methylomonas_MAG_cov

## Step 2.4 Calculate latesummer_day2Methylomonas_MAG_cov dict for each year    
os.mkdir("virus_n_MAG_tax_association/Methylomonas_virus_n_MAG")   
for year in year2IMGs:
    IMGs = year2IMGs[year]
    latesummer_day2Methylomonas_MAG_cov = {} # int(latesummer_day) => Methylomonas_MAG_cov
    for IMG in IMGs:
        latesummer_day = int(IMG2latesummer_day[IMG])
        Methylomonas_MAG_cov = IMG2Methylomonas_MAG_cov[IMG]
        latesummer_day2Methylomonas_MAG_cov[latesummer_day] = Methylomonas_MAG_cov

    ## Write down the result
    f = open(f"virus_n_MAG_tax_association/Methylomonas_virus_n_MAG/{year}.Methylomonas_MAG_cov.txt", 'w')
    f.write('latesummer_day\tMethylomonas_MAG_cov\n')
    for latesummer_day in sorted(latesummer_day2Methylomonas_MAG_cov.keys()):
        line = str(latesummer_day) + '\t' + str(latesummer_day2Methylomonas_MAG_cov[latesummer_day])
        f.write(line + '\n') 
    f.close()    
    
    
# Step 3 Get pmoC-containing Methylomonas virus abundance   
Species = {}  # $gn_rep => $gns
AMG_summary = {}  # $pro => $ko
KOs = {}  # $ko => 1
AMG_containing_viral_gn = {}  # $gn => 1
Old_gene2new_gene_map = {}  # $gene_old => $gene_new
PmoC_containing_viral_gn = {}  # $gn => 1

## Step 3.1 Store species info
with open("/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt", "r") as file:
    for line in file:
        line = line.strip()
        tmp = line.split("\t")
        gn_rep = tmp[0]
        gns = tmp[1]
        Species[gn_rep] = gns

## Step 3.2 Store AMG KO information
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

## Step 3.3 Get pmoC-containing viral genome
for pro in sorted(AMG_summary.keys()):
    gn = pro.split("__Ga")[0]
    ko = AMG_summary[pro]
    if gn in Species and ko == "K10946":
        PmoC_containing_viral_gn[gn] = 1
        
## Step 3.4 Store viral_gn2IMG2cov dict
viral_gn2IMG2cov = defaultdict(dict)
Header2 = []  # Store the header line
with open('MetaPop/Viral_gn2IMG2cov_norm_filtered.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Head'):
            Header2 = line.split('\t')
        else:
            tmp = line.split('\t')
            for i in range(1, len(tmp)):
                viral_gn = tmp[0]
                IMG = Header2[i]
                cov = tmp[i]
                if cov == 'NA':
                    cov = '0'
                viral_gn2IMG2cov[viral_gn][IMG] = float(cov)

## Step 3.5 Store viral_gn2host_tax dict and viral_gn2Methylomonas dict
viral_gn2host_tax = {} # viral_gn => host_tax
viral_gn2Methylomonas = set() # Set to store viral_gn with host as Methylomonas
with open('/storage1/data11/TYMEFLIES_phage/Host_prediction/Viral_gn2host_tax_final.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        viral_gn, host_tax = tmp[0], tmp[1]
        host_genus = host_tax.split(';s__', 1)[0]
        viral_gn2host_tax[viral_gn] = host_tax
        
        # Check if the host is Methylomonas and add the viral_gn to the set
        if host_genus in Methylomonas_genus_list:
            viral_gn2Methylomonas.add(viral_gn)    

## Step 3.6 Calculate IMG2pmoC_containing_Methylomonas_viral_gn_cov dict
IMG2pmoC_containing_Methylomonas_viral_gn_cov = {} # IMG => pmoC_containing_Methylomonas_viral_gn_cov (The sum cov of pmoC-containing Methylomonas viral gn (rep))            
for IMG in IMG2date:
    pmoC_containing_Methylomonas_viral_gn_cov = 0
    for viral_gn in viral_gn2IMG2cov:
        if viral_gn in PmoC_containing_viral_gn and viral_gn in viral_gn2Methylomonas:
            cov = viral_gn2IMG2cov[viral_gn][IMG]
            pmoC_containing_Methylomonas_viral_gn_cov += float(cov)
    IMG2pmoC_containing_Methylomonas_viral_gn_cov[IMG] = pmoC_containing_Methylomonas_viral_gn_cov    
                 
## Step 3.7 Calculate latesummer_day2pmoC_containing_Methylomonas_viral_gn_cov dicts for each year
for year in year2IMGs:
    IMGs = year2IMGs[year]
    latesummer_day2pmoC_containing_Methylomonas_viral_gn_cov = {} # int(latesummer_day) => pmoC_containing_Methylomonas_viral_gn_cov
    for IMG in IMGs:
        latesummer_day = int(IMG2latesummer_day[IMG])
        pmoC_containing_Methylomonas_viral_gn_cov = IMG2pmoC_containing_Methylomonas_viral_gn_cov[IMG]
        latesummer_day2pmoC_containing_Methylomonas_viral_gn_cov[latesummer_day] = pmoC_containing_Methylomonas_viral_gn_cov

    ## Write down the result
    f = open(f"virus_n_MAG_tax_association/Methylomonas_virus_n_MAG/{year}.pmoC_containing_Methylomonas_viral_gn_cov.txt", 'w')
    f.write('latesummer_day\tpmoC_containing_Methylomonas_viral_gn_cov\n')
    for latesummer_day in sorted(latesummer_day2pmoC_containing_Methylomonas_viral_gn_cov.keys()):
        line = str(latesummer_day) + '\t' + str(latesummer_day2pmoC_containing_Methylomonas_viral_gn_cov[latesummer_day])
        f.write(line + '\n')
    f.close()    


# Step 4 Get no-pmoC-containing Methylomonas virus abundance 
## Step 4.1 Store non-AMG containing virus to IMG 2 cov norm filtered
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
                
## Step 4.2 Store no-pmoC-containing Methylomonas viral gn set
no_pmoC_containing_Methylomonas_viral_gn_set = set()
for viral_gn in viral_gn2Methylomonas:
    if viral_gn in no_AMG_viral_gn2IMG2cov_norm_filtered:
        no_pmoC_containing_Methylomonas_viral_gn_set.add(viral_gn)
    elif viral_gn in viral_gn2IMG2cov and viral_gn not in PmoC_containing_viral_gn:
        no_pmoC_containing_Methylomonas_viral_gn_set.add(viral_gn)        
        
## Step 4.3 Calculate IMG2no_pmoC_containing_Methylomonas_viral_gn_cov dict
IMG2no_pmoC_containing_Methylomonas_viral_gn_cov = {} # IMG => no_pmoC_containing_Methylomonas_viral_gn_cov          
for IMG in IMG2date:
    no_pmoC_containing_Methylomonas_viral_gn_cov = 0
    for viral_gn in no_AMG_viral_gn2IMG2cov_norm_filtered:
        if viral_gn in no_pmoC_containing_Methylomonas_viral_gn_set:
            cov = no_AMG_viral_gn2IMG2cov_norm_filtered[viral_gn][IMG]
            no_pmoC_containing_Methylomonas_viral_gn_cov += cov
    for viral_gn in viral_gn2IMG2cov:
        if viral_gn in no_pmoC_containing_Methylomonas_viral_gn_set:
            cov = viral_gn2IMG2cov[viral_gn][IMG]
            no_pmoC_containing_Methylomonas_viral_gn_cov += cov
    IMG2no_pmoC_containing_Methylomonas_viral_gn_cov[IMG] = no_pmoC_containing_Methylomonas_viral_gn_cov 
     
## Step 4.4 Calculate latesummer_day2no_pmoC_containing_Methylomonas_viral_gn_cov dicts for each year
for year in year2IMGs:
    IMGs = year2IMGs[year]
    latesummer_day2no_pmoC_containing_Methylomonas_viral_gn_cov = {} # int(latesummer_day) => no_pmoC_containing_Methylomonas_viral_gn_cov
    for IMG in IMGs:
        latesummer_day = int(IMG2latesummer_day[IMG])
        no_pmoC_containing_Methylomonas_viral_gn_cov = IMG2no_pmoC_containing_Methylomonas_viral_gn_cov[IMG]
        latesummer_day2no_pmoC_containing_Methylomonas_viral_gn_cov[latesummer_day] = no_pmoC_containing_Methylomonas_viral_gn_cov

    ## Write down the result
    f = open(f"virus_n_MAG_tax_association/Methylomonas_virus_n_MAG/{year}.no_pmoC_containing_Methylomonas_viral_gn_cov.txt", 'w')
    f.write('latesummer_day\tno_pmoC_containing_Methylomonas_viral_gn_cov\n')
    for latesummer_day in sorted(latesummer_day2no_pmoC_containing_Methylomonas_viral_gn_cov.keys()):
        line = str(latesummer_day) + '\t' + str(latesummer_day2no_pmoC_containing_Methylomonas_viral_gn_cov[latesummer_day])
        f.write(line + '\n')  
    f.close()        