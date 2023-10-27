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
    
# Aim: (1) The percentage of species rep containing AMG  
#      (2) Calculate psbA-containing species rep to species with any members containing psbA ratio 
#      (3) Calculate virus completeness to specific AMG-containing species members AMG containing percentage

# Step 1 Get the species rep that contains AMG
## Step 1.1 Store species info
Species = {} # gn_rep => gns (for example, "3300033816__vRhyme_113,3300034101__vRhyme_unbinned4136,3300042863__vRhyme_unbinned973")
with open("/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt", "r") as file:
    for line in file:
        line = line.strip()
        tmp = line.split("\t")
        gn_rep = tmp[0]
        gns = tmp[1]
        Species[gn_rep] = gns

## Step 1.2 Store AMG KO information
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

## Step 1.3 Get the percentage of species rep containing AMG
## The stringent result 
species_rep_containing_AMG = []
for species_rep in Species:
    if species_rep in AMG_containing_viral_gn:
        species_rep_containing_AMG.append(species_rep)

percentage_species_rep_containing_AMG = 0
percentage_species_rep_containing_AMG = len(species_rep_containing_AMG) / len(Species.keys())  
print(percentage_species_rep_containing_AMG)    

## The loose result
species_rep_containing_AMG2 = set() # Containing the species rep if any members within the species contain AMG
for species_rep in Species:
    gn_list = Species[species_rep].split(',')
    for gn in gn_list:
        if gn in AMG_containing_viral_gn:
            species_rep_containing_AMG2.add(species_rep)

percentage_species_rep_containing_AMG2 = 0
percentage_species_rep_containing_AMG2 = len(species_rep_containing_AMG2) / len(Species.keys())  
print(percentage_species_rep_containing_AMG2) 
               

# Step 2 Calculate psbA-containing species rep to species with any members containing psbA ratio 
psbA_containing_species_rep = []
for species_rep in Species:
    if species_rep in AMG_containing_viral_gn:
        if 'K02703' in AMG_containing_viral_gn[species_rep]:
            psbA_containing_species_rep.append(species_rep)
            
species_with_any_members_containing_psbA = set()            
for species_rep in Species:
    gn_list = Species[species_rep].split(',')
    for gn in gn_list:
        if gn in AMG_containing_viral_gn and 'K02703' in AMG_containing_viral_gn[gn]:
            species_with_any_members_containing_psbA.add(gn)
                
ratio = len(psbA_containing_species_rep) / len(species_with_any_members_containing_psbA)
print(f"psbA-containing species rep to species with any members containing psbA ratio is {ratio}") 


# Step 3 Calculate virus completeness to specific AMG-containing species members AMG containing percentage
AMG2KO = {
    'psbA': 'K02703',
    'psbD': 'K02706',
    'pmoC': 'K10946',
    'ahbD': 'K22227',
    'katG': 'K03782',
    'gpx': 'K00432',
    'cobS': 'K09882',
    'cobT': 'K09883',
    'nadE': 'K01916',
    'cysC': 'K00860',
    'cysD': 'K00957',
    'cysH': 'K00390',
    'cysK': 'K01738'
}

viral_gn2completeness = defaultdict(list) # viral_gn => [0] checkv_quality [1] completeness [2] completeness_item
with open('/storage1/data11/TYMEFLIES_phage/CheckV_phage_bin_all/quality_summary.tsv', 'r') as lines:
    for line in lines:
        if not line.startswith('contig'):
            tmp = line.split('\t')
            viral_gn, checkv_quality, completeness = tmp[0], tmp[7], tmp[9]
            viral_gn2completeness[viral_gn].append(checkv_quality)
            viral_gn2completeness[viral_gn].append(completeness)
            if completeness == 'NA':
                viral_gn2completeness[viral_gn].append('NA') 
            elif float(completeness) <= 100 and float(completeness) > 75:
                viral_gn2completeness[viral_gn].append('75-100')           
            elif float(completeness) <= 75 and float(completeness) > 50:
                viral_gn2completeness[viral_gn].append('50-75')   
            elif float(completeness) <= 50 and float(completeness) > 25:
                viral_gn2completeness[viral_gn].append('25-50')   
            elif float(completeness) <= 25 and float(completeness) >= 0:
                viral_gn2completeness[viral_gn].append('0-25')
                
checkv_quality_category = ['Complete', 'High-quality', 'Medium-quality', 'Low-quality', 'Not-determined']
completeness_item_category = ['75-100', '50-75', '25-50', '0-25', 'NA']

for AMG in AMG2KO:
    print(f"{AMG}:{AMG2KO[AMG]}")
    ## Get the species rep containing specific AMG and set the minimal species size to be n
    n = 1 # The minimal species size
    species_rep_containing_specific_AMG = set()
    for species_rep in Species:
        gn_list = Species[species_rep].split(',')
        if len(gn_list) >= n:
            if species_rep in AMG_containing_viral_gn and AMG2KO[AMG] in AMG_containing_viral_gn[species_rep]:
                species_rep_containing_specific_AMG.add(species_rep)
                
    ## Get the whole viruses from species_rep_containing_specific_AMG list
    virus_list = [] # Store the whole virus list
    for species_rep in species_rep_containing_specific_AMG:
        gn_list = Species[species_rep].split(',')
        for item in gn_list:
            virus_list.append(item) 

    ## Calculate virus completeness to AMG containing percentage          
    completeness_item_category2AMG_containing_percentage = {} # completeness_item => AMG_containing_percentage
    for completeness_item in completeness_item_category:
        AMG_containing_percentage = 0
        viral_gn_of_this_category = []
        for viral_gn in virus_list:
            if viral_gn2completeness[viral_gn][2] == completeness_item:
                viral_gn_of_this_category.append(viral_gn)
                
        viral_gn_of_this_category_containing_AMG = []
        for viral_gn in viral_gn_of_this_category:
            if viral_gn in AMG_containing_viral_gn and AMG2KO[AMG] in AMG_containing_viral_gn[viral_gn]:
                viral_gn_of_this_category_containing_AMG.append(viral_gn) 
      
        if len(viral_gn_of_this_category) > 0:
            AMG_containing_percentage = len(viral_gn_of_this_category_containing_AMG) / len(viral_gn_of_this_category)
            completeness_item_category2AMG_containing_percentage[completeness_item] = str(AMG_containing_percentage) + "|" + str(len(viral_gn_of_this_category))
    print(completeness_item_category2AMG_containing_percentage)
        