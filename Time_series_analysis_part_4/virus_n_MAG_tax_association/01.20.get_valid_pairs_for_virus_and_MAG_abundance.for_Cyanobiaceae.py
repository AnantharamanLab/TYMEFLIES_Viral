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
    
# Aim: Generate the input table for plotting box plots for virus and host abundance ratios of three groups
# The order of abundances are:
# For Cyanobiaceae:
# (1) Cyanobiaceae_MAG (2) psbA_containing_Cyanobiaceae_viral_gn (3) no_psbA_containing_Cyanobiaceae_viral_gn


# Step 1 Store Cyanobiaceae_MAG abundance from interpolation results
## Step 1.1 Store Cyanobiaceae_MAG_year2highest_y_value dict
Cyanobiaceae_MAG_year2highest_y_value = {} # year => highest_y_value
with open('Cyanobiaceae_virus_n_MAG/Cyanobiaceae_MAG.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        Cyanobiaceae_MAG_year2highest_y_value[year] = highest_y_value

## Step 1.2 Store Cyanobiaceae_MAG_date2y_value dict
Cyanobiaceae_MAG_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
Cyanobiaceae_MAG_interpolation_result_files = glob('Cyanobiaceae_virus_n_MAG/*.Cyanobiaceae_MAG.interpolation_results.txt')
for Cyanobiaceae_MAG_interpolation_result_file in Cyanobiaceae_MAG_interpolation_result_files:
    year = Cyanobiaceae_MAG_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(Cyanobiaceae_MAG_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            Cyanobiaceae_MAG_date2y_value[f"{year}|{date}"] = y_value
            
  
# Step 2 Store psbA_containing_Cyanobiaceae_viral_gn abundance from interpolation results
## Step 2.1 Store psbA_containing_Cyanobiaceae_viral_gn_year2highest_y_value dict
psbA_containing_Cyanobiaceae_viral_gn_year2highest_y_value = {} # year => highest_y_value
with open('Cyanobiaceae_virus_n_MAG/psbA_containing_Cyanobiaceae_viral_gn.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        psbA_containing_Cyanobiaceae_viral_gn_year2highest_y_value[year] = highest_y_value

## Step 2.2 Store psbA_containing_Cyanobiaceae_viral_gn_date2y_value dict
psbA_containing_Cyanobiaceae_viral_gn_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
psbA_containing_Cyanobiaceae_viral_gn_interpolation_result_files = glob('Cyanobiaceae_virus_n_MAG/*.psbA_containing_Cyanobiaceae_viral_gn.interpolation_results.txt')
for psbA_containing_Cyanobiaceae_viral_gn_interpolation_result_file in psbA_containing_Cyanobiaceae_viral_gn_interpolation_result_files:
    year = psbA_containing_Cyanobiaceae_viral_gn_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(psbA_containing_Cyanobiaceae_viral_gn_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            psbA_containing_Cyanobiaceae_viral_gn_date2y_value[f"{year}|{date}"] = y_value
            

# Step 3 Store no_psbA_containing_Cyanobiaceae_viral_gn abundance from interpolation results
## Step 3.1 Store no_psbA_containing_Cyanobiaceae_viral_gn_year2highest_y_value dict
no_psbA_containing_Cyanobiaceae_viral_gn_year2highest_y_value = {} # year => highest_y_value
with open('Cyanobiaceae_virus_n_MAG/no_psbA_containing_Cyanobiaceae_viral_gn.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        no_psbA_containing_Cyanobiaceae_viral_gn_year2highest_y_value[year] = highest_y_value

## Step 3.2 Store no_psbA_containing_Cyanobiaceae_viral_gn_date2y_value dict
no_psbA_containing_Cyanobiaceae_viral_gn_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
no_psbA_containing_Cyanobiaceae_viral_gn_interpolation_result_files = glob('Cyanobiaceae_virus_n_MAG/*.no_psbA_containing_Cyanobiaceae_viral_gn.interpolation_results.txt')
for no_psbA_containing_Cyanobiaceae_viral_gn_interpolation_result_file in no_psbA_containing_Cyanobiaceae_viral_gn_interpolation_result_files:
    year = no_psbA_containing_Cyanobiaceae_viral_gn_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(no_psbA_containing_Cyanobiaceae_viral_gn_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            no_psbA_containing_Cyanobiaceae_viral_gn_date2y_value[f"{year}|{date}"] = y_value

            
# Step 4 Filter valid pairs
date2valid_pair = {} # date => [Cyanobiaceae_MAG_abun, psbA_containing_Cyanobiaceae_viral_gn_abun, no_psbA_containing_Cyanobiaceae_viral_gn_abun]
for date in Cyanobiaceae_MAG_date2y_value:
    Cyanobiaceae_MAG_abun_ratio = Cyanobiaceae_MAG_date2y_value[date]
    psbA_containing_Cyanobiaceae_viral_gn_abun_ratio = psbA_containing_Cyanobiaceae_viral_gn_date2y_value[date]
    no_psbA_containing_Cyanobiaceae_viral_gn_abun_ratio = no_psbA_containing_Cyanobiaceae_viral_gn_date2y_value[date]
    
    year = date.split('|')[0]    
    if Cyanobiaceae_MAG_abun_ratio != 'nan' and psbA_containing_Cyanobiaceae_viral_gn_abun_ratio != 'nan' and no_psbA_containing_Cyanobiaceae_viral_gn_abun_ratio != 'nan':
        if Cyanobiaceae_MAG_year2highest_y_value[year] != 'nan' and psbA_containing_Cyanobiaceae_viral_gn_year2highest_y_value[year] != 'nan' and no_psbA_containing_Cyanobiaceae_viral_gn_year2highest_y_value[year] != 'nan': 
            Cyanobiaceae_MAG_abun = float(Cyanobiaceae_MAG_abun_ratio) * float(Cyanobiaceae_MAG_year2highest_y_value[year]) / 100    
            psbA_containing_Cyanobiaceae_viral_gn_abun = float(psbA_containing_Cyanobiaceae_viral_gn_abun_ratio) * float(psbA_containing_Cyanobiaceae_viral_gn_year2highest_y_value[year]) / 100 
            no_psbA_containing_Cyanobiaceae_viral_gn_abun = float(no_psbA_containing_Cyanobiaceae_viral_gn_abun_ratio) * float(no_psbA_containing_Cyanobiaceae_viral_gn_year2highest_y_value[year]) / 100                     
            if float(Cyanobiaceae_MAG_abun_ratio) > 10 and float(psbA_containing_Cyanobiaceae_viral_gn_abun_ratio) > 10 and float(no_psbA_containing_Cyanobiaceae_viral_gn_abun_ratio) > 10:
                date2valid_pair[date] = [str(Cyanobiaceae_MAG_abun), str(psbA_containing_Cyanobiaceae_viral_gn_abun), str(no_psbA_containing_Cyanobiaceae_viral_gn_abun)]
            
f = open('Cyanobiaceae_virus_n_MAG/psbA_date2valid_pair.for_Cyanobiaceae.txt', 'w') 
for date in date2valid_pair:
    valid_pair = date2valid_pair[date]
    line = date + '\t' + '\t'.join(valid_pair) + '\n'
    f.write(line)
f.close()    
    
           
            