#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    import datetime
    from collections import defaultdict
    import seaborn as sns
    import pandas as pd
    from tabulate import tabulate
    import matplotlib.pyplot as plt
    import matplotlib as mpl
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: Get seasonal KO-carrying viral genome distribution


# Step 1 Store AMG summary information
ko2ko_detail = {} # ko => ko_detail 
ko2IMG_set = defaultdict(set) # ko => set(IMG); This set contains the IMGs that this KO-carrying viral genome can be found within
with open('AMG_analysis/AMG_summary.txt', 'r') as file:
    for line in file:
        line = line.strip()
        if not line.startswith("Pro"):
            tmp = line.split("\t")
            pro = tmp[0]
            ko = tmp[2]
            ko_detail = tmp[3]
            IMG = pro.split('__')[0]
            ko2ko_detail[ko] = ko_detail
            ko2IMG_set[ko].add(IMG)
file.close()


# Step 2 Store season2IMG_no and IMG2season dict
season2IMG_no = defaultdict(int) # season => IMG_no; The number of IMGs that belong to the season
IMG2season = {} # IMG => season
with open('TYMEFLIES_metagenome_info.txt', 'r') as lines:
    for line in lines:
        line = line.strip('\n')
        if not line.startswith("IMG"):
            tmp = line.split("\t")
            season = tmp[10]
            season2IMG_no[season] += 1
            IMG = tmp[0]
            IMG2season[IMG] = season
lines.close()            


# Step 3 Store High Occurrence KO list
high_occurrence_ko_list = []
with open('AMG_analysis/KO2occurrence_n_abundance.txt', 'r') as lines:
    for line in lines:
        line = line.strip('\n')
        tmp = line.split('\t')
        if int(tmp[1]) >= 400:
            high_occurrence_ko_list.append(tmp[0])
lines.close() 


# Step 4 Get the ko2season2ko_carrying_genome_distribution_percentage dict
ko2season2ko_carrying_genome_distribution_percentage = defaultdict(dict) # ko => season => ko_carrying_genome_distribution_percentage
## Note: Only high occurrence ko was included
for ko in high_occurrence_ko_list:
    IMG_set_from_ko = ko2IMG_set[ko] # The IMG set from ko
    for season in season2IMG_no:
        IMG_no = season2IMG_no[season] # The number of IMGs that belong to the season 
        IMG_no_containing_ko = 0 # The number of IMGs that contain the ko
        for IMG in IMG_set_from_ko:
            if IMG2season[IMG] == season:
                IMG_no_containing_ko += 1
        ko_carrying_genome_distribution_percentage = IMG_no_containing_ko / IMG_no     
        ko2season2ko_carrying_genome_distribution_percentage[ko][season] = ko_carrying_genome_distribution_percentage


# Step 5 Write down the result
season_list = ["Spring", "Clearwater", "Early Summer", "Late Summer", "Fall", "Ice-on"]
with open('AMG_analysis/KO2season2KO_carrying_genome_distribution_percentage.txt', "w") as file:
    ## Write table header
    header = "\t".join(season_list)
    file.write("KO\t" + header + "\n")

    ## Write table rows
    for ko in high_occurrence_ko_list:
        row_values = [str(ko2season2ko_carrying_genome_distribution_percentage.get(ko, {}).get(season, '')) for season in season_list]
        row = ko + "\t" + "\t".join(row_values) + "\n"
        file.write(row)
        