#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    from datetime import datetime  
    from collections import defaultdict
    import statistics
    import subprocess
    import scipy.stats
    import statsmodels.stats.multitest as multi
    import numpy as np  
    from pathlib import Path  
    from scipy.interpolate import interp1d    
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
# Aim: Conduct correlation analysis between each pair of virus/host and env parameters
# Note: Use the new method - try to find whether there are correlations between the number  
#       of days maintained at > 20% of the peak abundance (both cyanobacteria and cyanoviruses) 
#       and the mean environmental parameters during the summer seasons (both Early Summer and Late Summer)


def calc_over_20perc_days(perc_list, date_list):
    days = 0 # Store the day numbers that have over 20 percent bacteria/viral abundance
    perc_chunks_list = []  # Initialize an empty list to store the chunks of percentages
    date_chunks_list = []  # Initialize an empty list to store the chunks of dates

    perc_chunk = []  # Initialize an empty sublist to store the current perc_chunk
    date_chunk = []  # Initialize an empty sublist to store the current date_chunk

    for i, item in enumerate(perc_list):
        if item != 'nan' and float(item) > 20:
            perc_chunk.append(item)
            date_chunk.append(date_list[i])
        else:
            if perc_chunk:  # Check if the current chunk is not empty
                perc_chunks_list.append(perc_chunk)  # Append the current chunk to the result
                date_chunks_list.append(date_chunk)  # Append the corresponding chunk
                perc_chunk = []  # Reset the current perc_chunk
                date_chunk = []  # Reset the current date_chunk

    # Check if the last chunk is not empty and add it to the result
    if perc_chunk:
        perc_chunks_list.append(perc_chunk)
        date_chunks_list.append(date_chunk)

    for date_chunk in date_chunks_list:
        days += int(date_chunk[-1]) - int(date_chunk[0])

    return days


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
    
    
# Step 2 Get the number of days for psbA-containing Cyanobacteria virus maintained at > 20% of the peak abundance
Species = {}  # $gn_rep => $gns
AMG_summary = {}  # $pro => $ko
KOs = {}  # $ko => 1
AMG_containing_viral_gn = {}  # $gn => 1
Old_gene2new_gene_map = {}  # $gene_old => $gene_new
PsbA_containing_viral_gn = {}  # $gn => 1

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

## Step 2.3 Get psbA-containing viral genome
for pro in sorted(AMG_summary.keys()):
    gn = pro.split("__Ga")[0]
    ko = AMG_summary[pro]
    if gn in Species and ko == "K02703":
        PsbA_containing_viral_gn[gn] = 1
              
## Step 2.4 Store viral_gn2IMG2cov dict
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

# Step 2.5 Store viral_gn2host_tax dict and viral_gn2Cyanobacteria dict
viral_gn2host_tax = {} # viral_gn => host_tax
viral_gn2Cyanobacteria = set() # Set to store viral_gn with host as Cyanobacteria
with open('/storage1/data11/TYMEFLIES_phage/Host_prediction/Viral_gn2host_tax_final.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        viral_gn, host_tax = tmp[0], tmp[1]
        host_phylum = host_tax.split(';c__', 1)[0]
        viral_gn2host_tax[viral_gn] = host_tax
        
        # Check if the host is Cyanobacteria and add the viral_gn to the set
        if host_phylum == 'd__Bacteria;p__Cyanobacteria':
            viral_gn2Cyanobacteria.add(viral_gn)  

## Step 2.6 Calculate IMG2psbA_containing_Cyanobacteria_viral_gn_cov dict
IMG2psbA_containing_Cyanobacteria_viral_gn_cov = {} # IMG => psbA_containing_Cyanobacteria_viral_gn_cov (The sum cov of psbA-containing Cyanobacteria viral gn (rep))            
for IMG in IMG2date:
    psbA_containing_Cyanobacteria_viral_gn_cov = 0
    for viral_gn in viral_gn2IMG2cov:
        if viral_gn in PsbA_containing_viral_gn and viral_gn in viral_gn2Cyanobacteria:
            cov = viral_gn2IMG2cov[viral_gn][IMG]
            psbA_containing_Cyanobacteria_viral_gn_cov += cov
    IMG2psbA_containing_Cyanobacteria_viral_gn_cov[IMG] = psbA_containing_Cyanobacteria_viral_gn_cov    
                 
## Step 2.7 Calculate earlysummer_day2psbA_containing_Cyanobacteria_viral_gn_cov dicts for each year
for year in year2IMGs:
    IMGs = year2IMGs[year]
    earlysummer_day2psbA_containing_Cyanobacteria_viral_gn_cov = {} # int(earlysummer_day) => psbA_containing_Cyanobacteria_viral_gn_cov
    for IMG in IMGs:
        earlysummer_day = int(IMG2earlysummer_day[IMG])
        psbA_containing_Cyanobacteria_viral_gn_cov = IMG2psbA_containing_Cyanobacteria_viral_gn_cov[IMG]
        earlysummer_day2psbA_containing_Cyanobacteria_viral_gn_cov[earlysummer_day] = psbA_containing_Cyanobacteria_viral_gn_cov

    ### Write down the result
    f = open(f"virus_n_host_n_env_association/{year}.psbA_containing_Cyanobacteria_viral_gn_cov.txt", 'w')
    f.write('earlysummer_day\tpsbA_containing_Cyanobacteria_viral_gn_cov\n')
    for earlysummer_day in sorted(earlysummer_day2psbA_containing_Cyanobacteria_viral_gn_cov.keys()):
        line = str(earlysummer_day) + '\t' + str(earlysummer_day2psbA_containing_Cyanobacteria_viral_gn_cov[earlysummer_day])
        f.write(line + '\n')
    f.close()   

## Step 2.8 Store the psbA-containing Cyanobacteria viral gn cov files for each year
psbA_containing_Cyanobacteria_viral_gn_cov_files = glob("virus_n_host_n_env_association/*.psbA_containing_Cyanobacteria_viral_gn_cov.txt")
year2psbA_containing_Cyanobacteria_viral_gn_cov = defaultdict(dict)
for psbA_containing_Cyanobacteria_viral_gn_cov_file in psbA_containing_Cyanobacteria_viral_gn_cov_files:
    year = Path(psbA_containing_Cyanobacteria_viral_gn_cov_file).stem.split('.')[0]    
    x_list = []
    y_list = []
    with open(psbA_containing_Cyanobacteria_viral_gn_cov_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('earlysummer'):
                tmp = line.split('\t')
                x_list.append(float(tmp[0]))  # Convert x value to float
                y_list.append(float(tmp[1]))  # Convert y value to float
    year2psbA_containing_Cyanobacteria_viral_gn_cov[year]['x'] = x_list 
    year2psbA_containing_Cyanobacteria_viral_gn_cov[year]['y'] = y_list       

## Step 2.9 Create interpolation functions for each year's dataset
interpolation_functions = {}
for year, data in year2psbA_containing_Cyanobacteria_viral_gn_cov.items():
    x = np.array(data['x'])  # Convert x values to numpy array
    y = np.array(data['y'])  # Convert y values to numpy array
    interpolation_functions[year] = interp1d(x, y, bounds_error=False)

### Define mean line x points
x_mean = np.arange(-45, 160, 5)

### Evaluate interpolation functions at x_mean for each year
interpolation_results = defaultdict(dict)
for year, interp_func in interpolation_functions.items():
    y_mean = interp_func(x_mean)
    interpolation_results[year]['x_mean'] = x_mean
    interpolation_results[year]['y_mean'] = y_mean

### Store the interpolation results in text files
output_directory = "virus_n_host_n_env_association"  # Replace with your desired output directory

year2highest_y_value = {} # year => highest_y_value
for year, data in interpolation_results.items():
    output_filename = f"{year}.psbA_containing_Cyanobacteria_viral_gn.interpolation_results.txt"
    output_file_path = Path(output_directory) / output_filename

    x_mean_values = data['x_mean']
    y_mean_values = data['y_mean']

    # Find the highest value in y_mean_values
    highest_y_value = max(y_mean_values)
    year2highest_y_value[year] = highest_y_value
    
    with open(output_file_path, 'w') as output_file:
        for x_val, y_val in zip(x_mean_values, y_mean_values):
            # Calculate the percentage of y_val
            percentage_y_val = (y_val / highest_y_value) * 100

            # Print the modified y_val
            output_file.write(f"{x_val}\t{percentage_y_val:.2f}\n")  

## Step 2.10 Get the number of days for psbA-containing Cyanobacteria virus maintained at > 20% of the peak abundance for each year
year2psbA_containing_Cyanobacteria_viral_gn_abun_maintained_over_20perc_days = {} # year => over_20perc_days
psbA_containing_Cyanobacteria_viral_gn_interpolation_results_files = glob("virus_n_host_n_env_association/*.psbA_containing_Cyanobacteria_viral_gn.interpolation_results.txt")
for file in psbA_containing_Cyanobacteria_viral_gn_interpolation_results_files:
    year = Path(file).stem.split('.', 1)[0]
    perc_list = []
    date_list = []
    with open(file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            date, perc = line.split('\t')
            perc_list.append(perc)
            date_list.append(date)
            
    over_20perc_days = calc_over_20perc_days(perc_list, date_list)
    year2psbA_containing_Cyanobacteria_viral_gn_abun_maintained_over_20perc_days[year] = over_20perc_days
   
## Step 2.11 Write down the year2psbA_containing_Cyanobacteria_viral_gn_abun_maintained_over_20perc_days dict
f = open('virus_n_host_n_env_association/year2psbA_containing_Cyanobacteria_viral_gn_abun_maintained_over_20perc_days.txt', 'w')
for year in sorted(year2psbA_containing_Cyanobacteria_viral_gn_abun_maintained_over_20perc_days.keys()):
    days = year2psbA_containing_Cyanobacteria_viral_gn_abun_maintained_over_20perc_days[year]
    line = year + '\t' + str(days)
    f.write(line + '\n')
 

# Step 3 Get the number of days for no-psbA-containing Cyanobacteria virus maintained at > 20% of the peak abundance 
## Step 3.1 Store no-AMG containing virus to IMG 2 cov norm filtered
no_AMG_viral_gn2IMG2cov_norm_filtered = defaultdict(dict) # viral_gn => IMG => cov_norm_filtered
Header2 = [] # Store the header line
with open('MetaPop/no_AMG_viral_gn2IMG2cov_norm_filtered.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Head'):
            tmp = line.split('\t')
            Header2 = tmp
        else:
            tmp = line.split('\t')
            viral_gn = tmp[0]
            for i in range(1, len(tmp)):
                IMG = Header2[i]
                cov_norm_filtered = tmp[i]
                no_AMG_viral_gn2IMG2cov_norm_filtered[viral_gn][IMG] = float(cov_norm_filtered)   
     
## Step 3.2 Store no-psbA-containing Cyanobacteria viral gn set
no_psbA_containing_Cyanobacteria_viral_gn_set = set()
for viral_gn in viral_gn2Cyanobacteria:
    if viral_gn in no_AMG_viral_gn2IMG2cov_norm_filtered:
        no_psbA_containing_Cyanobacteria_viral_gn_set.add(viral_gn)
    elif viral_gn in viral_gn2IMG2cov and viral_gn not in PsbA_containing_viral_gn:
        no_psbA_containing_Cyanobacteria_viral_gn_set.add(viral_gn)

## Step 3.3 Calculate IMG2no_psbA_containing_Cyanobacteria_viral_gn_cov dict
IMG2no_psbA_containing_Cyanobacteria_viral_gn_cov = {} # IMG => no_psbA_containing_Cyanobacteria_viral_gn_cov          
for IMG in IMG2date:
    no_psbA_containing_Cyanobacteria_viral_gn_cov = 0
    for viral_gn in no_AMG_viral_gn2IMG2cov_norm_filtered:
        if viral_gn in no_psbA_containing_Cyanobacteria_viral_gn_set:
            cov = no_AMG_viral_gn2IMG2cov_norm_filtered[viral_gn][IMG]
            no_psbA_containing_Cyanobacteria_viral_gn_cov += cov
    for viral_gn in viral_gn2IMG2cov:
        if viral_gn in no_psbA_containing_Cyanobacteria_viral_gn_set:
            cov = viral_gn2IMG2cov[viral_gn][IMG]
            no_psbA_containing_Cyanobacteria_viral_gn_cov += cov            
    IMG2no_psbA_containing_Cyanobacteria_viral_gn_cov[IMG] = no_psbA_containing_Cyanobacteria_viral_gn_cov
     
## Step 3.4 Calculate earlysummer_day2no_psbA_containing_Cyanobacteria_viral_gn_cov dicts for each year
for year in year2IMGs:
    IMGs = year2IMGs[year]
    earlysummer_day2no_psbA_containing_Cyanobacteria_viral_gn_cov = {} # int(earlysummer_day) => no_psbA_containing_Cyanobacteria_viral_gn_cov
    for IMG in IMGs:
        earlysummer_day = int(IMG2earlysummer_day[IMG])
        no_psbA_containing_Cyanobacteria_viral_gn_cov = IMG2no_psbA_containing_Cyanobacteria_viral_gn_cov[IMG]
        earlysummer_day2no_psbA_containing_Cyanobacteria_viral_gn_cov[earlysummer_day] = no_psbA_containing_Cyanobacteria_viral_gn_cov

    ### Write down the result
    f = open(f"virus_n_host_n_env_association/{year}.no_psbA_containing_Cyanobacteria_viral_gn_cov.txt", 'w')
    f.write('earlysummer_day\tno_psbA_containing_Cyanobacteria_viral_gn_cov\n')
    for earlysummer_day in sorted(earlysummer_day2no_psbA_containing_Cyanobacteria_viral_gn_cov.keys()):
        line = str(earlysummer_day) + '\t' + str(earlysummer_day2no_psbA_containing_Cyanobacteria_viral_gn_cov[earlysummer_day])
        f.write(line + '\n')  
    f.close() 
    
## Step 3.5 Store the no-psbA-containing Cyanobacteria viral gn cov files for each year
no_psbA_containing_Cyanobacteria_viral_gn_cov_files = glob("virus_n_host_n_env_association/*.no_psbA_containing_Cyanobacteria_viral_gn_cov.txt")
year2no_psbA_containing_Cyanobacteria_viral_gn_cov = defaultdict(dict)
for no_psbA_containing_Cyanobacteria_viral_gn_cov_file in no_psbA_containing_Cyanobacteria_viral_gn_cov_files:
    year = Path(no_psbA_containing_Cyanobacteria_viral_gn_cov_file).stem.split('.')[0]    
    x_list = []
    y_list = []
    with open(no_psbA_containing_Cyanobacteria_viral_gn_cov_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('earlysummer'):
                tmp = line.split('\t')
                x_list.append(float(tmp[0]))  # Convert x value to float
                y_list.append(float(tmp[1]))  # Convert y value to float
    year2no_psbA_containing_Cyanobacteria_viral_gn_cov[year]['x'] = x_list 
    year2no_psbA_containing_Cyanobacteria_viral_gn_cov[year]['y'] = y_list       

## Step 3.6 Create interpolation functions for each year's dataset
interpolation_functions = {}
for year, data in year2no_psbA_containing_Cyanobacteria_viral_gn_cov.items():
    x = np.array(data['x'])  # Convert x values to numpy array
    y = np.array(data['y'])  # Convert y values to numpy array
    interpolation_functions[year] = interp1d(x, y, bounds_error=False)

### Define mean line x points
x_mean = np.arange(-45, 160, 5)

### Evaluate interpolation functions at x_mean for each year
interpolation_results = defaultdict(dict)
for year, interp_func in interpolation_functions.items():
    y_mean = interp_func(x_mean)
    interpolation_results[year]['x_mean'] = x_mean
    interpolation_results[year]['y_mean'] = y_mean

### Store the interpolation results in text files
output_directory = "virus_n_host_n_env_association"  # Replace with your desired output directory

year2highest_y_value = {} # year => highest_y_value
for year, data in interpolation_results.items():
    output_filename = f"{year}.no_psbA_containing_Cyanobacteria_viral_gn.interpolation_results.txt"
    output_file_path = Path(output_directory) / output_filename

    x_mean_values = data['x_mean']
    y_mean_values = data['y_mean']

    # Find the highest value in y_mean_values
    highest_y_value = max(y_mean_values)
    year2highest_y_value[year] = highest_y_value
    
    with open(output_file_path, 'w') as output_file:
        for x_val, y_val in zip(x_mean_values, y_mean_values):
            # Calculate the percentage of y_val
            percentage_y_val = (y_val / highest_y_value) * 100

            # Print the modified y_val
            output_file.write(f"{x_val}\t{percentage_y_val:.2f}\n")  

## Step 3.7 Get the number of days for no-psbA-containing Cyanobacteria virus maintained at > 20% of the peak abundance for each year
year2no_psbA_containing_Cyanobacteria_viral_gn_abun_maintained_over_20perc_days = {} # year => over_20perc_days
no_psbA_containing_Cyanobacteria_viral_gn_interpolation_results_files = glob("virus_n_host_n_env_association/*.no_psbA_containing_Cyanobacteria_viral_gn.interpolation_results.txt")
for file in no_psbA_containing_Cyanobacteria_viral_gn_interpolation_results_files:
    year = Path(file).stem.split('.', 1)[0]
    perc_list = []
    date_list = []
    with open(file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            date, perc = line.split('\t')
            perc_list.append(perc)
            date_list.append(date)
            
    over_20perc_days = calc_over_20perc_days(perc_list, date_list)
    year2no_psbA_containing_Cyanobacteria_viral_gn_abun_maintained_over_20perc_days[year] = over_20perc_days
   
## Step 3.8 Write down the year2no_psbA_containing_Cyanobacteria_viral_gn_abun_maintained_over_20perc_days dict
f = open('virus_n_host_n_env_association/year2no_psbA_containing_Cyanobacteria_viral_gn_abun_maintained_over_20perc_days.txt', 'w')
for year in sorted(year2no_psbA_containing_Cyanobacteria_viral_gn_abun_maintained_over_20perc_days.keys()):
    days = year2no_psbA_containing_Cyanobacteria_viral_gn_abun_maintained_over_20perc_days[year]
    line = year + '\t' + str(days)
    f.write(line + '\n')  


# Step 4 Get the number of days for all Cyanobacteria virus (including both psbA-containing and no-psbA-containing viruses) maintained at > 20% of the peak abundance  
## Step 4.1 Store the Cyanobacteria viral gn cov files for each year
no_psbA_containing_Cyanobacteria_viral_gn_cov_files = glob("virus_n_host_n_env_association/*.no_psbA_containing_Cyanobacteria_viral_gn_cov.txt")
psbA_containing_Cyanobacteria_viral_gn_cov_files = glob("virus_n_host_n_env_association/*.psbA_containing_Cyanobacteria_viral_gn_cov.txt")
year2Cyanobacteria_viral_gn_cov = defaultdict(dict)
for no_psbA_containing_Cyanobacteria_viral_gn_cov_file in no_psbA_containing_Cyanobacteria_viral_gn_cov_files:
    year = Path(no_psbA_containing_Cyanobacteria_viral_gn_cov_file).stem.split('.')[0]    
    x_list1 = []
    y_list1 = []
    x_list2 = []
    y_list2 = []    
    psbA_containing_Cyanobacteria_viral_gn_cov_file = no_psbA_containing_Cyanobacteria_viral_gn_cov_file.replace('no_psbA', 'psbA', 1)
    with open(no_psbA_containing_Cyanobacteria_viral_gn_cov_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('earlysummer'):
                tmp = line.split('\t')
                x_list1.append(float(tmp[0]))  # Convert x value to float
                y_list1.append(float(tmp[1]))  # Convert y value to float
    with open(psbA_containing_Cyanobacteria_viral_gn_cov_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('earlysummer'):
                tmp = line.split('\t')
                x_list2.append(float(tmp[0]))  # Convert x value to float
                y_list2.append(float(tmp[1]))  # Convert y value to float 
    
    y_list = [] # The combined list of two y_lists
    # Check if both lists have the same length
    if len(y_list1) == len(y_list2):
        y_list = [a + b for a, b in zip(y_list1, y_list2)]
    else:
        print("Lists have different lengths")                
    year2Cyanobacteria_viral_gn_cov[year]['x'] = x_list1
    year2Cyanobacteria_viral_gn_cov[year]['y'] = y_list       

## Step 4.2 Create interpolation functions for each year's dataset
interpolation_functions = {}
for year, data in year2Cyanobacteria_viral_gn_cov.items():
    x = np.array(data['x'])  # Convert x values to numpy array
    y = np.array(data['y'])  # Convert y values to numpy array
    interpolation_functions[year] = interp1d(x, y, bounds_error=False)

### Define mean line x points
x_mean = np.arange(-45, 160, 5)

### Evaluate interpolation functions at x_mean for each year
interpolation_results = defaultdict(dict)
for year, interp_func in interpolation_functions.items():
    y_mean = interp_func(x_mean)
    interpolation_results[year]['x_mean'] = x_mean
    interpolation_results[year]['y_mean'] = y_mean

### Store the interpolation results in text files
output_directory = "virus_n_host_n_env_association"  # Replace with your desired output directory

year2highest_y_value = {} # year => highest_y_value
for year, data in interpolation_results.items():
    output_filename = f"{year}.Cyanobacteria_viral_gn.interpolation_results.txt"
    output_file_path = Path(output_directory) / output_filename

    x_mean_values = data['x_mean']
    y_mean_values = data['y_mean']

    # Find the highest value in y_mean_values
    highest_y_value = max(y_mean_values)
    year2highest_y_value[year] = highest_y_value
    
    with open(output_file_path, 'w') as output_file:
        for x_val, y_val in zip(x_mean_values, y_mean_values):
            # Calculate the percentage of y_val
            percentage_y_val = (y_val / highest_y_value) * 100

            # Print the modified y_val
            output_file.write(f"{x_val}\t{percentage_y_val:.2f}\n")  

## Step 4.3 Get the number of days for all Cyanobacteria virus maintained at > 20% of the peak abundance for each year
year2Cyanobacteria_viral_gn_abun_maintained_over_20perc_days = {} # year => over_20perc_days
Cyanobacteria_viral_gn_interpolation_results_files = glob("virus_n_host_n_env_association/*.Cyanobacteria_viral_gn.interpolation_results.txt")
for file in Cyanobacteria_viral_gn_interpolation_results_files:
    year = Path(file).stem.split('.', 1)[0]
    perc_list = []
    date_list = []
    with open(file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            date, perc = line.split('\t')
            perc_list.append(perc)
            date_list.append(date)
            
    over_20perc_days = calc_over_20perc_days(perc_list, date_list)
    year2Cyanobacteria_viral_gn_abun_maintained_over_20perc_days[year] = over_20perc_days
   
## Step 4.4 Write down the year2Cyanobacteria_viral_gn_abun_maintained_over_20perc_days dict
f = open('virus_n_host_n_env_association/year2Cyanobacteria_viral_gn_abun_maintained_over_20perc_days.txt', 'w')
for year in sorted(year2Cyanobacteria_viral_gn_abun_maintained_over_20perc_days.keys()):
    days = year2Cyanobacteria_viral_gn_abun_maintained_over_20perc_days[year]
    line = year + '\t' + str(days)
    f.write(line + '\n')  

    
# Step 5 Get the number of days for Cyanobacteria maintained at > 20% of the peak abundance     
## Step 5.1 Store MAG2IMG2abun dict
MAG2IMG2abun = defaultdict(dict)
Header3 = []  # Store the header line
with open('virus_n_MAG_tax_association/MAG2IMG2abun.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('head'):
            Header3 = line.split('\t')
        else:
            tmp = line.split('\t')
            for i in range(1, len(tmp)):
                MAG = tmp[0]
                IMG = Header3[i]
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
            if 'p__Cyanobacteria' in tax:
                Cyanobacteria_MAG.append(MAG)              
        
## Step 5.3 Calculate IMG2Cyanobacteria_MAG_cov dict          
IMG2Cyanobacteria_MAG_cov = {} # IMG => Cyanobacteria_MAG_cov (The sum cov of Cyanobacteria MAGs)            
for IMG in IMG2date:
    Cyanobacteria_MAG_cov = 0
    for MAG in Cyanobacteria_MAG:
        cov = MAG2IMG2abun[MAG][IMG]
        Cyanobacteria_MAG_cov += cov
    IMG2Cyanobacteria_MAG_cov[IMG] = Cyanobacteria_MAG_cov

## Step 5.4 Calculate earlysummer_day2Cyanobacteria_MAG_cov dict for each year 
for year in year2IMGs:
    IMGs = year2IMGs[year]
    earlysummer_day2Cyanobacteria_MAG_cov = {} # int(earlysummer_day) => Cyanobacteria_MAG_cov
    for IMG in IMGs:
        earlysummer_day = int(IMG2earlysummer_day[IMG])
        Cyanobacteria_MAG_cov = IMG2Cyanobacteria_MAG_cov[IMG]
        earlysummer_day2Cyanobacteria_MAG_cov[earlysummer_day] = Cyanobacteria_MAG_cov

    ### Write down the result
    f = open(f"virus_n_host_n_env_association/{year}.Cyanobacteria_MAG_cov.txt", 'w')
    f.write('earlysummer_day\tCyanobacteria_MAG_cov\n')
    for earlysummer_day in sorted(earlysummer_day2Cyanobacteria_MAG_cov.keys()):
        line = str(earlysummer_day) + '\t' + str(earlysummer_day2Cyanobacteria_MAG_cov[earlysummer_day])
        f.write(line + '\n') 
    f.close() 
    
## Step 5.5 Store the Cyanobacteria MAG cov files for each year
Cyanobacteria_MAG_cov_files = glob("virus_n_host_n_env_association/*.Cyanobacteria_MAG_cov.txt")
year2Cyanobacteria_MAG_cov = defaultdict(dict)
for Cyanobacteria_MAG_cov_file in Cyanobacteria_MAG_cov_files:
    year = Path(Cyanobacteria_MAG_cov_file).stem.split('.')[0]    
    x_list = []
    y_list = []
    with open(Cyanobacteria_MAG_cov_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('earlysummer'):
                tmp = line.split('\t')
                x_list.append(float(tmp[0]))  # Convert x value to float
                y_list.append(float(tmp[1]))  # Convert y value to float
    year2Cyanobacteria_MAG_cov[year]['x'] = x_list 
    year2Cyanobacteria_MAG_cov[year]['y'] = y_list       

## Step 5.6 Create interpolation functions for each year's dataset
interpolation_functions = {}
for year, data in year2Cyanobacteria_MAG_cov.items():
    x = np.array(data['x'])  # Convert x values to numpy array
    y = np.array(data['y'])  # Convert y values to numpy array
    interpolation_functions[year] = interp1d(x, y, bounds_error=False)

### Define mean line x points
x_mean = np.arange(-45, 160, 5)

### Evaluate interpolation functions at x_mean for each year
interpolation_results = defaultdict(dict)
for year, interp_func in interpolation_functions.items():
    y_mean = interp_func(x_mean)
    interpolation_results[year]['x_mean'] = x_mean
    interpolation_results[year]['y_mean'] = y_mean

### Store the interpolation results in text files
output_directory = "virus_n_host_n_env_association"  # Replace with your desired output directory

year2highest_y_value = {} # year => highest_y_value
for year, data in interpolation_results.items():
    output_filename = f"{year}.Cyanobacteria_MAG.interpolation_results.txt"
    output_file_path = Path(output_directory) / output_filename

    x_mean_values = data['x_mean']
    y_mean_values = data['y_mean']

    # Find the highest value in y_mean_values
    highest_y_value = max(y_mean_values)
    year2highest_y_value[year] = highest_y_value
    
    with open(output_file_path, 'w') as output_file:
        for x_val, y_val in zip(x_mean_values, y_mean_values):
            # Calculate the percentage of y_val
            percentage_y_val = (y_val / highest_y_value) * 100

            # Print the modified y_val
            output_file.write(f"{x_val}\t{percentage_y_val:.2f}\n")   
                        
## Step 5.7 Get the number of days for Cyanobacteria maintained at > 20% of the peak abundance for each year
year2Cyanobacteria_MAG_abun_maintained_over_20perc_days = {} # year => over_20perc_days
Cyanobacteria_MAG_interpolation_results_files = glob("virus_n_host_n_env_association/*.Cyanobacteria_MAG.interpolation_results.txt")
for file in Cyanobacteria_MAG_interpolation_results_files:
    year = Path(file).stem.split('.', 1)[0]
    perc_list = []
    date_list = []
    with open(file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            date, perc = line.split('\t')
            perc_list.append(perc)
            date_list.append(date)
            
    over_20perc_days = calc_over_20perc_days(perc_list, date_list)
    year2Cyanobacteria_MAG_abun_maintained_over_20perc_days[year] = over_20perc_days
   
## Step 5.8 Write down the year2Cyanobacteria_MAG_abun_maintained_over_20perc_days dict
f = open('virus_n_host_n_env_association/year2Cyanobacteria_MAG_abun_maintained_over_20perc_days.txt', 'w')
for year in sorted(year2Cyanobacteria_MAG_abun_maintained_over_20perc_days.keys()):
    days = year2Cyanobacteria_MAG_abun_maintained_over_20perc_days[year]
    line = year + '\t' + str(days)
    f.write(line + '\n')             
    