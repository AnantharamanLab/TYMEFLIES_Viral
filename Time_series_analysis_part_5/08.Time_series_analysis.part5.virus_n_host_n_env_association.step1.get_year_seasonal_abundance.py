#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    import datetime
    from collections import defaultdict
    import statistics
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: Get the year-season abundance for:
#      (1) psbA AMG cov (parse intermediate results from "AMG_gene2IMG2cov_ratio.txt" and "Viral_gn2IMG2cov_norm_filtered.txt")
#      (2) psbA-containing Cyanobacteria virus abundance
#      (3) no-psbA-containing Cyanobacteria virus abundance  
#      (4) Cyanobacteria abundance


# Step 1 Store AMG KO information and year-season information
## Step 1.1 Store AMG KO information
AMG_summary = {}  # pro => ko
KOs = {}  # ko => ko_detail
IMG2date = {}  # IMG => date_n_season
AMG_containing_viral_gn = {}  # gn => 1
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

            IMG = pro.split('__', 1)[0]
            IMG2date[IMG] = date_n_season
            KOs[ko] = ko_detail

            gn = pro.rsplit('__', 1)[0]
            AMG_containing_viral_gn[gn] = 1

### Change the old gene to new gene
Old_gene2new_gene_map = {}  # gene_old => gene_new
with open("New_gene2old_gene_map.txt", "r") as file:
    for line in file:
        line = line.strip()
        tmp = line.split("\t")
        gene_new = tmp[0]
        gene_old = tmp[1]
        Old_gene2new_gene_map[gene_old] = gene_new

keys_to_delete = []

### Create a list of keys to delete
for pro in AMG_summary:
    if pro in Old_gene2new_gene_map:
        keys_to_delete.append(pro)

### Iterate over the list of keys to delete and make the changes
for pro in keys_to_delete:
    ko = AMG_summary[pro]
    gene_new = Old_gene2new_gene_map[pro]
    del AMG_summary[pro]  # Delete the old gene and its value
    AMG_summary[gene_new] = ko  # Add the new gene and its value

## Step 1.2 Store season metagenome information
Season2num = {}  # season => num_metagenome
Season2IMG = defaultdict(list)  # season => [IMGs]
Year_season2num = {}  # year_season => num_metagenome
Year_season2IMG = defaultdict(list)  # year_season => [IMGs]
with open("TYMEFLIES_metagenome_info.txt", "r") as file:
    for line in file:
        line = line.strip()
        if not line.startswith("IMG"):
            tmp = line.split("\t")
            IMG = tmp[0]
            date = tmp[8]
            season = tmp[10]
            Season2num[season] = Season2num.get(season, 0) + 1

            Season2IMG[season].append(IMG)    

            year = tmp[8].split("-")[0]
            year_season = f"{year}-{season}"
            Year_season2num[year_season] = Year_season2num.get(year_season, 0) + 1

            Year_season2IMG[year_season].append(IMG)    

Season = ['Spring', 'Clearwater', 'Early Summer', 'Late Summer', 'Fall', 'Ice-on']
Year = [str(year) for year in range(2000, 2020)]
Year_season = [f"{year}-{season}" for year in Year for season in Season]


# Step 2 Get psbA_AMG2year_season2cov
## Step 2.1 Store "AMG_gene2IMG2cov_ratio.txt"
AMG_gene2IMG2cov_ratio = defaultdict(dict)  # AMG => IMG => cov_ratio
Header = []  # Store the header line
IMG_set = set()
with open('MetaPop/AMG_gene2IMG2cov_ratio.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith("Head"):
            Header = line.split("\t")
        else:
            tmp = line.split("\t")
            AMG = tmp[0]
            for i in range(1, len(tmp)):
                IMG = Header[i]
                IMG_set.add(IMG)
                cov_ratio = tmp[i]
                if cov_ratio == 'Both viral gn and AMG absent' or cov_ratio == 'Viral gn absent, AMG present':
                    cov_ratio = 'NA'
                AMG_gene2IMG2cov_ratio[AMG][IMG] = cov_ratio 

## Step 2.2 Store "Viral_gn2IMG2cov_norm_filtered.txt"
Viral_gn2IMG2cov_norm_filtered = defaultdict(dict)  # viral_gn => IMG => cov_norm_filtered
Header2 = []  # Store the header line
with open('MetaPop/Viral_gn2IMG2cov_norm_filtered.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith("Head"):
            Header2 = line.split("\t")
        else:
            tmp = line.split("\t")
            viral_gn = tmp[0]
            for i in range(1, len(tmp)):
                IMG = Header2[i]
                cov_norm_filtered = tmp[i]
                Viral_gn2IMG2cov_norm_filtered[viral_gn][IMG] = cov_norm_filtered

## Step 2.3 Get AMG_gene2IMG2cov dict
AMG_gene2IMG2cov = defaultdict(dict)
for AMG in AMG_gene2IMG2cov_ratio:
    for IMG in IMG_set:
        cov_ratio = AMG_gene2IMG2cov_ratio[AMG][IMG] # cov_ratio of AMG to viral_gn
        viral_gn = AMG.rsplit('__', 1)[0]
        viral_gn_cov_norm_filtered = Viral_gn2IMG2cov_norm_filtered[viral_gn][IMG] # the cov_norm_filtered of the corresponding viral_gn
        
        cov = 'NA'
        if cov_ratio != 'NA' and viral_gn_cov_norm_filtered != 'NA':
            cov = float(cov_ratio) * float(viral_gn_cov_norm_filtered)
        AMG_gene2IMG2cov[AMG][IMG] = cov
        
## Step 2.4 Get IMG2all_psbA_AMG_cov dict
IMG2all_psbA_AMG_cov = {} # IMG => all_psbA_AMG_cov
for IMG in IMG_set:
    cov = 0 # The cov value for all_psbA_AMG in this IMG
    for AMG in AMG_gene2IMG2cov_ratio:        
        if AMG_summary[AMG] == 'K02703': # Check if the AMG is a psbA_AMG            
            if AMG_gene2IMG2cov[AMG][IMG] != 'NA':
                cov += float(AMG_gene2IMG2cov[AMG][IMG])
    IMG2all_psbA_AMG_cov[IMG] = cov  
       
## Step 2.5 Get year_season2all_psbA_AMG_cov
year_season2all_psbA_AMG_cov = {} # year_season => all_psbA_AMG_cov 
for year_season in Year_season:
    IMG_list = Year_season2IMG[year_season] # Store all metagenomes from this year_season
    Cov_collection = [] # Store all values for all_psbA_AMG_cov
    cov_for_this_year_season = 0  # Store the mean all_psbA_AMG_cov for this year_season
    for IMG in IMG_list:
        cov = IMG2all_psbA_AMG_cov[IMG]
        Cov_collection.append(cov)   

    if len(Cov_collection) > 0:
        cov_for_this_year_season = sum(Cov_collection) / len(Cov_collection)
                
    year_season2all_psbA_AMG_cov[year_season] = cov_for_this_year_season
                
## Step 2.6 Write down year_season2all_psbA_AMG_cov
### Header line
header = "Head\t" + "\t".join(Year_season)

### Initialize the table content
table_content = [header]

### Write down the result line
line = 'all_psbA_AMG_cov'
for year_season in Year_season:
    cov = year_season2all_psbA_AMG_cov[year_season]
    line += '\t' + str(cov)
table_content.append(line)    

### Write the table to the file
#os.mkdir("virus_n_host_n_env_association")
with open("virus_n_host_n_env_association/year_season2all_psbA_AMG_cov.txt", "w") as outfile:
    outfile.write("\n".join(table_content))   


# Step 3 Get psbA-containing Cyanobacteria virus abundance
## Step 3.1 Store species info
Species = {}  # $gn_rep => $gns
with open("/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt", "r") as file:
    for line in file:
        line = line.strip()
        tmp = line.split("\t")
        gn_rep = tmp[0]
        gns = tmp[1]
        Species[gn_rep] = gns
        
## Step 3.2 Get psbA-containing viral genome
PsbA_containing_viral_gn = {} # gn => 1; Store the psbA-containg viral gn (species rep)
for pro in sorted(AMG_summary.keys()):
    gn = pro.split("__Ga")[0]
    ko = AMG_summary[pro]
    if gn in Species and ko == "K02703":
        PsbA_containing_viral_gn[gn] = 1

## Step 3.3 Store viral_gn2IMG2cov dict
viral_gn2IMG2cov = defaultdict(dict) # Here, viral_gn contain all AMG-containing viral species rep gn
Header3 = []  # Store the header line
with open('MetaPop/Viral_gn2IMG2cov_norm_filtered.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Head'):
            Header3 = line.split('\t')
        else:
            tmp = line.split('\t')
            for i in range(1, len(tmp)):
                viral_gn = tmp[0]
                IMG = Header3[i]
                cov = tmp[i]
                viral_gn2IMG2cov[viral_gn][IMG] = cov # Here, cov is a string, it can be 'NA' or number
                
## Step 3.4 Store viral_gn2host_tax dict and viral_gn2Cyanobacteria dict
viral_gn2host_tax = {} # viral_gn => host_tax
viral_gn2Cyanobacteria = set() # Set to store viral_gn with host as Cyanobacteria
with open('/storage1/data11/TYMEFLIES_phage/Host_prediction/Viral_gn2host_tax_final.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        viral_gn, host_tax = tmp[0], tmp[1]
        viral_gn2host_tax[viral_gn] = host_tax
        
        # Check if the host_tax is Cyanobacteria and add the viral_gn to the set
        if 'd__Bacteria;p__Cyanobacteria' in host_tax:
            viral_gn2Cyanobacteria.add(viral_gn)
            
## Step 3.5 Get IMG2all_psbA_containing_Cyanobacteria_viral_gn_abun dict
IMG2all_psbA_containing_Cyanobacteria_viral_gn_abun = {} # IMG => all_psbA_containing_Cyanobacteria_viral_gn_abun
for IMG in IMG_set:
    abun = 0 # The abun for all_psbA_containing_Cyanobacteria_viral_gn
    for viral_gn in PsbA_containing_viral_gn: # This viral gn (species rep) has psbA 
        if viral_gn in viral_gn2Cyanobacteria: # if this viral_gn (species rep) is Cyanobacteria virus
            if viral_gn2IMG2cov[viral_gn][IMG] != 'NA':
                abun += float(viral_gn2IMG2cov[viral_gn][IMG])
    IMG2all_psbA_containing_Cyanobacteria_viral_gn_abun[IMG] = abun

## Step 3.6 Parse to get year_season2all_psbA_containing_Cyanobacteria_viral_gn_abun dict
year_season2all_psbA_containing_Cyanobacteria_viral_gn_abun = {} # year_season => all_psbA_containing_Cyanobacteria_viral_gn_abun 
for year_season in Year_season:
    IMG_list = Year_season2IMG[year_season] # Store all metagenomes from this year_season
    Abun_collection = [] # Store all values for all_psbA_containing_Cyanobacteria_viral_gn_abun
    abun_for_this_year_season = 0  # Store the mean all_psbA_containing_Cyanobacteria_viral_gn_abun for this year_season
    for IMG in IMG_list:
        abun = IMG2all_psbA_containing_Cyanobacteria_viral_gn_abun[IMG]
        Abun_collection.append(abun)   

    if len(Abun_collection) > 0:
        abun_for_this_year_season = sum(Abun_collection) / len(Abun_collection)
                
    year_season2all_psbA_containing_Cyanobacteria_viral_gn_abun[year_season] = abun_for_this_year_season

## Step 3.7 Write down year_season2all_psbA_containing_Cyanobacteria_viral_gn_abun
### Header line
header2 = "Head\t" + "\t".join(Year_season)

### Initialize the table content
table_content2 = [header2]

### Write down the result line
line2 = 'all_psbA_containing_Cyanobacteria_viral_gn_abun'
for year_season in Year_season:
    abun = year_season2all_psbA_containing_Cyanobacteria_viral_gn_abun[year_season]
    line2 += '\t' + str(abun)
table_content2.append(line2)    

### Write the table to the file
#os.mkdir("virus_n_host_n_env_association")
with open("virus_n_host_n_env_association/year_season2all_psbA_containing_Cyanobacteria_viral_gn_abun.txt", "w") as outfile:
    outfile.write("\n".join(table_content2))      


# Step 4 Get no-psbA-containing Cyanobacteria virus abundance
## Step 4.1 Store no-AMG containing virus to IMG 2 cov norm filtered
no_AMG_viral_gn2IMG2cov_norm_filtered = defaultdict(dict) # viral_gn => IMG => cov_norm_filtered
Header4 = [] # Store the header line
with open('MetaPop/no_AMG_viral_gn2IMG2cov_norm_filtered.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Head'):
            tmp = line.split('\t')
            Header4 = tmp
        else:
            tmp = line.split('\t')
            viral_gn = tmp[0]
            for i in range(1, len(tmp)):
                IMG = Header4[i]
                cov_norm_filtered = tmp[i]
                no_AMG_viral_gn2IMG2cov_norm_filtered[viral_gn][IMG] = float(cov_norm_filtered)  

## Step 4.2 Store no-psbA-containing Cyanobacteria viral gn set
no_psbA_containing_Cyanobacteria_viral_gn_set = set()
for viral_gn in viral_gn2Cyanobacteria:
    if viral_gn in no_AMG_viral_gn2IMG2cov_norm_filtered:
        no_psbA_containing_Cyanobacteria_viral_gn_set.add(viral_gn)
    elif viral_gn in viral_gn2IMG2cov and viral_gn not in PsbA_containing_viral_gn:
        no_psbA_containing_Cyanobacteria_viral_gn_set.add(viral_gn)
        
## Step 4.3 Get IMG2all_no_psbA_containing_Cyanobacteria_viral_gn_abun dict
IMG2all_no_psbA_containing_Cyanobacteria_viral_gn_abun = {} # IMG => no_psbA_containing_Cyanobacteria_viral_gn_abun
for IMG in IMG_set:
    abun = 0 # The abun for all_no_psbA_containing_Cyanobacteria_viral_gn
    for viral_gn in no_psbA_containing_Cyanobacteria_viral_gn_set:
        if viral_gn in viral_gn2IMG2cov:
            if viral_gn2IMG2cov[viral_gn][IMG] != 'NA':
                abun += float(viral_gn2IMG2cov[viral_gn][IMG])
        elif viral_gn in no_AMG_viral_gn2IMG2cov_norm_filtered:
            abun += float(no_AMG_viral_gn2IMG2cov_norm_filtered[viral_gn][IMG])              
    IMG2all_no_psbA_containing_Cyanobacteria_viral_gn_abun[IMG] = abun        

## Step 4.4 Parse to get year_season2all_no_psbA_containing_Cyanobacteria_viral_gn_abun dict
year_season2all_no_psbA_containing_Cyanobacteria_viral_gn_abun = {} # year_season => all_no_psbA_containing_Cyanobacteria_viral_gn_abun 
for year_season in Year_season:
    IMG_list = Year_season2IMG[year_season] # Store all metagenomes from this year_season
    Abun_collection = [] # Store all values for all_no_psbA_containing_Cyanobacteria_viral_gn_abun
    abun_for_this_year_season = 0  # Store the mean all_no_psbA_containing_Cyanobacteria_viral_gn_abun for this year_season
    for IMG in IMG_list:
        abun = IMG2all_no_psbA_containing_Cyanobacteria_viral_gn_abun[IMG]
        Abun_collection.append(abun)   

    if len(Abun_collection) > 0:
        abun_for_this_year_season = sum(Abun_collection) / len(Abun_collection)
                
    year_season2all_no_psbA_containing_Cyanobacteria_viral_gn_abun[year_season] = abun_for_this_year_season   

## Step 4.5 Write down year_season2all_no_psbA_containing_Cyanobacteria_viral_gn_abun
### Header line
header3 = "Head\t" + "\t".join(Year_season)

### Initialize the table content
table_content3 = [header3]

### Write down the result line
line3 = 'all_no_psbA_containing_Cyanobacteria_viral_gn_abun'
for year_season in Year_season:
    abun = year_season2all_no_psbA_containing_Cyanobacteria_viral_gn_abun[year_season]
    line3 += '\t' + str(abun)
table_content3.append(line3)    

### Write the table to the file
#os.mkdir("virus_n_host_n_env_association")
with open("virus_n_host_n_env_association/year_season2all_no_psbA_containing_Cyanobacteria_viral_gn_abun.txt", "w") as outfile:
    outfile.write("\n".join(table_content3))         


# Step 5 Get Cyanobacteria abundance
## Step 5.1 Store MAG2IMG2abun dict
MAG2IMG2abun = defaultdict(dict)
Header5 = []  # Store the header line
with open('virus_n_MAG_tax_association/MAG2IMG2abun.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('head'):
            Header5 = line.split('\t')
        else:
            tmp = line.split('\t')
            for i in range(1, len(tmp)):
                MAG = tmp[0]
                IMG = Header5[i]
                abun = tmp[i]
                MAG2IMG2abun[MAG][IMG] = float(abun)
                
## Step 5.2 Store Cyanobacteria_MAG list
Cyanobacteria_MAG = []
with open('/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.rep_MAG.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('IMG'):
            tmp = line.split('\t')
            MAG, tax = tmp[0], tmp[1]
            if 'd__Bacteria;p__Cyanobacteria' in tax:
                Cyanobacteria_MAG.append(MAG) 

## Step 5.3 Get IMG2all_Cyanobacteria_MAG_abun dict
IMG2all_Cyanobacteria_MAG_abun = {} 
for IMG in IMG_set:
    abun = 0 # The abun for all_Cyanobacteria_MAG
    for MAG in Cyanobacteria_MAG:
        abun += float(MAG2IMG2abun[MAG][IMG])
    IMG2all_Cyanobacteria_MAG_abun[IMG] = abun    
                
## Step 5.4 Parse to get year_season2all_Cyanobacteria_MAG_abun dict
year_season2all_Cyanobacteria_MAG_abun = {} # year_season => all_Cyanobacteria_MAG_abun 
for year_season in Year_season:
    IMG_list = Year_season2IMG[year_season] # Store all metagenomes from this year_season
    Abun_collection = [] # Store all values for all_Cyanobacteria_MAG_abun
    abun_for_this_year_season = 0  # Store the mean all_Cyanobacteria_MAG_abun for this year_season
    for IMG in IMG_list:
        abun = IMG2all_Cyanobacteria_MAG_abun[IMG]
        Abun_collection.append(abun)   

    if len(Abun_collection) > 0:
        abun_for_this_year_season = sum(Abun_collection) / len(Abun_collection)
                
    year_season2all_Cyanobacteria_MAG_abun[year_season] = abun_for_this_year_season    

## Step 5.5 Write down year_season2all_Cyanobacteria_MAG_abun
### Header line
header4 = "Head\t" + "\t".join(Year_season)

### Initialize the table content
table_content4 = [header4]

### Write down the result line
line4 = 'all_Cyanobacteria_MAG_abun'
for year_season in Year_season:
    abun = year_season2all_Cyanobacteria_MAG_abun[year_season]
    line4 += '\t' + str(abun)
table_content4.append(line4)    

### Write the table to the file
#os.mkdir("virus_n_host_n_env_association")
with open("virus_n_host_n_env_association/year_season2all_Cyanobacteria_MAG_abun.txt", "w") as outfile:
    outfile.write("\n".join(table_content4))       