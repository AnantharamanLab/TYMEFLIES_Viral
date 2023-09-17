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
# For Methylomonas:
# (1) Methylomonas_MAG (2) pmoC_containing_Methylomonas_viral_gn (3) no_pmoC_containing_Methylomonas_viral_gn


# Step 1 Store Methylomonas_MAG abundance from interpolation results
## Step 1.1 Store Methylomonas_MAG_year2highest_y_value dict
Methylomonas_MAG_year2highest_y_value = {} # year => highest_y_value
with open('Methylomonas_virus_n_MAG/Methylomonas_MAG.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        Methylomonas_MAG_year2highest_y_value[year] = highest_y_value

## Step 1.2 Store Methylomonas_MAG_date2y_value dict
Methylomonas_MAG_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
Methylomonas_MAG_interpolation_result_files = glob('Methylomonas_virus_n_MAG/*.Methylomonas_MAG.interpolation_results.txt')
for Methylomonas_MAG_interpolation_result_file in Methylomonas_MAG_interpolation_result_files:
    year = Methylomonas_MAG_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(Methylomonas_MAG_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            Methylomonas_MAG_date2y_value[f"{year}|{date}"] = y_value
            
  
# Step 2 Store pmoC_containing_Methylomonas_viral_gn abundance from interpolation results
## Step 2.1 Store pmoC_containing_Methylomonas_viral_gn_year2highest_y_value dict
pmoC_containing_Methylomonas_viral_gn_year2highest_y_value = {} # year => highest_y_value
with open('Methylomonas_virus_n_MAG/pmoC_containing_Methylomonas_viral_gn.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        pmoC_containing_Methylomonas_viral_gn_year2highest_y_value[year] = highest_y_value

## Step 2.2 Store pmoC_containing_Methylomonas_viral_gn_date2y_value dict
pmoC_containing_Methylomonas_viral_gn_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
pmoC_containing_Methylomonas_viral_gn_interpolation_result_files = glob('Methylomonas_virus_n_MAG/*.pmoC_containing_Methylomonas_viral_gn.interpolation_results.txt')
for pmoC_containing_Methylomonas_viral_gn_interpolation_result_file in pmoC_containing_Methylomonas_viral_gn_interpolation_result_files:
    year = pmoC_containing_Methylomonas_viral_gn_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(pmoC_containing_Methylomonas_viral_gn_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            pmoC_containing_Methylomonas_viral_gn_date2y_value[f"{year}|{date}"] = y_value
            

# Step 3 Store no_pmoC_containing_Methylomonas_viral_gn abundance from interpolation results
## Step 3.1 Store no_pmoC_containing_Methylomonas_viral_gn_year2highest_y_value dict
no_pmoC_containing_Methylomonas_viral_gn_year2highest_y_value = {} # year => highest_y_value
with open('Methylomonas_virus_n_MAG/no_pmoC_containing_Methylomonas_viral_gn.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        no_pmoC_containing_Methylomonas_viral_gn_year2highest_y_value[year] = highest_y_value

## Step 3.2 Store no_pmoC_containing_Methylomonas_viral_gn_date2y_value dict
no_pmoC_containing_Methylomonas_viral_gn_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
no_pmoC_containing_Methylomonas_viral_gn_interpolation_result_files = glob('Methylomonas_virus_n_MAG/*.no_pmoC_containing_Methylomonas_viral_gn.interpolation_results.txt')
for no_pmoC_containing_Methylomonas_viral_gn_interpolation_result_file in no_pmoC_containing_Methylomonas_viral_gn_interpolation_result_files:
    year = no_pmoC_containing_Methylomonas_viral_gn_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(no_pmoC_containing_Methylomonas_viral_gn_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            no_pmoC_containing_Methylomonas_viral_gn_date2y_value[f"{year}|{date}"] = y_value

            
# Step 4 Filter valid pairs
date2valid_pair = {} # date => [Methylomonas_MAG_abun, pmoC_containing_Methylomonas_viral_gn_abun, no_pmoC_containing_Methylomonas_viral_gn_abun]
for date in Methylomonas_MAG_date2y_value:
    Methylomonas_MAG_abun_ratio = Methylomonas_MAG_date2y_value[date]
    pmoC_containing_Methylomonas_viral_gn_abun_ratio = pmoC_containing_Methylomonas_viral_gn_date2y_value[date]
    no_pmoC_containing_Methylomonas_viral_gn_abun_ratio = no_pmoC_containing_Methylomonas_viral_gn_date2y_value[date]
    
    year = date.split('|')[0]    
    if Methylomonas_MAG_abun_ratio != 'nan' and pmoC_containing_Methylomonas_viral_gn_abun_ratio != 'nan' and no_pmoC_containing_Methylomonas_viral_gn_abun_ratio != 'nan':
        if Methylomonas_MAG_year2highest_y_value[year] != 'nan' and pmoC_containing_Methylomonas_viral_gn_year2highest_y_value[year] != 'nan' and no_pmoC_containing_Methylomonas_viral_gn_year2highest_y_value[year] != 'nan': 
            Methylomonas_MAG_abun = float(Methylomonas_MAG_abun_ratio) * float(Methylomonas_MAG_year2highest_y_value[year]) / 100    
            pmoC_containing_Methylomonas_viral_gn_abun = float(pmoC_containing_Methylomonas_viral_gn_abun_ratio) * float(pmoC_containing_Methylomonas_viral_gn_year2highest_y_value[year]) / 100 
            no_pmoC_containing_Methylomonas_viral_gn_abun = float(no_pmoC_containing_Methylomonas_viral_gn_abun_ratio) * float(no_pmoC_containing_Methylomonas_viral_gn_year2highest_y_value[year]) / 100                     
            if float(Methylomonas_MAG_abun_ratio) > 10 and float(pmoC_containing_Methylomonas_viral_gn_abun_ratio) > 10 and float(no_pmoC_containing_Methylomonas_viral_gn_abun_ratio) > 10:
                date2valid_pair[date] = [str(Methylomonas_MAG_abun), str(pmoC_containing_Methylomonas_viral_gn_abun), str(no_pmoC_containing_Methylomonas_viral_gn_abun)]
            
f = open('Methylomonas_virus_n_MAG/pmoC_date2valid_pair.for_Methylomonas.txt', 'w') 
for date in date2valid_pair:
    valid_pair = date2valid_pair[date]
    line = date + '\t' + '\t'.join(valid_pair) + '\n'
    f.write(line)
f.close()    
    
           
            