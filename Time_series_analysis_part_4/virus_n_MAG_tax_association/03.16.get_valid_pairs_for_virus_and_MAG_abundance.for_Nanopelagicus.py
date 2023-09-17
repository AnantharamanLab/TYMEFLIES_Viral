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
# For Nanopelagicus:
# (1) Nanopelagicus_MAG (2) katG_containing_Nanopelagicus_viral_gn (3) no_katG_containing_Nanopelagicus_viral_gn


# Step 1 Store Nanopelagicus_MAG abundance from interpolation results
## Step 1.1 Store Nanopelagicus_MAG_year2highest_y_value dict
Nanopelagicus_MAG_year2highest_y_value = {} # year => highest_y_value
with open('Nanopelagicus_virus_n_MAG/Nanopelagicus_MAG.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        Nanopelagicus_MAG_year2highest_y_value[year] = highest_y_value

## Step 1.2 Store Nanopelagicus_MAG_date2y_value dict
Nanopelagicus_MAG_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
Nanopelagicus_MAG_interpolation_result_files = glob('Nanopelagicus_virus_n_MAG/*.Nanopelagicus_MAG.interpolation_results.txt')
for Nanopelagicus_MAG_interpolation_result_file in Nanopelagicus_MAG_interpolation_result_files:
    year = Nanopelagicus_MAG_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(Nanopelagicus_MAG_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            Nanopelagicus_MAG_date2y_value[f"{year}|{date}"] = y_value
            
  
# Step 2 Store katG_containing_Nanopelagicus_viral_gn abundance from interpolation results
## Step 2.1 Store katG_containing_Nanopelagicus_viral_gn_year2highest_y_value dict
katG_containing_Nanopelagicus_viral_gn_year2highest_y_value = {} # year => highest_y_value
with open('Nanopelagicus_virus_n_MAG/katG_containing_Nanopelagicus_viral_gn.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        katG_containing_Nanopelagicus_viral_gn_year2highest_y_value[year] = highest_y_value

## Step 2.2 Store katG_containing_Nanopelagicus_viral_gn_date2y_value dict
katG_containing_Nanopelagicus_viral_gn_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
katG_containing_Nanopelagicus_viral_gn_interpolation_result_files = glob('Nanopelagicus_virus_n_MAG/*.katG_containing_Nanopelagicus_viral_gn.interpolation_results.txt')
for katG_containing_Nanopelagicus_viral_gn_interpolation_result_file in katG_containing_Nanopelagicus_viral_gn_interpolation_result_files:
    year = katG_containing_Nanopelagicus_viral_gn_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(katG_containing_Nanopelagicus_viral_gn_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            katG_containing_Nanopelagicus_viral_gn_date2y_value[f"{year}|{date}"] = y_value
            

# Step 3 Store no_katG_containing_Nanopelagicus_viral_gn abundance from interpolation results
## Step 3.1 Store no_katG_containing_Nanopelagicus_viral_gn_year2highest_y_value dict
no_katG_containing_Nanopelagicus_viral_gn_year2highest_y_value = {} # year => highest_y_value
with open('Nanopelagicus_virus_n_MAG/no_katG_containing_Nanopelagicus_viral_gn.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        no_katG_containing_Nanopelagicus_viral_gn_year2highest_y_value[year] = highest_y_value

## Step 3.2 Store no_katG_containing_Nanopelagicus_viral_gn_date2y_value dict
no_katG_containing_Nanopelagicus_viral_gn_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
no_katG_containing_Nanopelagicus_viral_gn_interpolation_result_files = glob('Nanopelagicus_virus_n_MAG/*.no_katG_containing_Nanopelagicus_viral_gn.interpolation_results.txt')
for no_katG_containing_Nanopelagicus_viral_gn_interpolation_result_file in no_katG_containing_Nanopelagicus_viral_gn_interpolation_result_files:
    year = no_katG_containing_Nanopelagicus_viral_gn_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(no_katG_containing_Nanopelagicus_viral_gn_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            no_katG_containing_Nanopelagicus_viral_gn_date2y_value[f"{year}|{date}"] = y_value

            
# Step 4 Filter valid pairs
date2valid_pair = {} # date => [Nanopelagicus_MAG_abun, no_katG_containing_Nanopelagicus_viral_gn_abun]
# Here, we only consider Nanopelagicus_MAG_abun and no_katG_containing_Nanopelagicus_viral_gn_abun for this case
for date in Nanopelagicus_MAG_date2y_value:
    Nanopelagicus_MAG_abun_ratio = Nanopelagicus_MAG_date2y_value[date]
    no_katG_containing_Nanopelagicus_viral_gn_abun_ratio = no_katG_containing_Nanopelagicus_viral_gn_date2y_value[date]
    
    year = date.split('|')[0]    
    if Nanopelagicus_MAG_abun_ratio != 'nan' and no_katG_containing_Nanopelagicus_viral_gn_abun_ratio != 'nan':
        if Nanopelagicus_MAG_year2highest_y_value[year] != 'nan' and no_katG_containing_Nanopelagicus_viral_gn_year2highest_y_value[year] != 'nan': 
            Nanopelagicus_MAG_abun = float(Nanopelagicus_MAG_abun_ratio) * float(Nanopelagicus_MAG_year2highest_y_value[year]) / 100    
            no_katG_containing_Nanopelagicus_viral_gn_abun = float(no_katG_containing_Nanopelagicus_viral_gn_abun_ratio) * float(no_katG_containing_Nanopelagicus_viral_gn_year2highest_y_value[year]) / 100                     
            if float(Nanopelagicus_MAG_abun_ratio) > 10 and float(no_katG_containing_Nanopelagicus_viral_gn_abun_ratio) > 10:
                date2valid_pair[date] = [str(Nanopelagicus_MAG_abun), str(no_katG_containing_Nanopelagicus_viral_gn_abun)]
            
f = open('Nanopelagicus_virus_n_MAG/katG_date2valid_pair.for_Nanopelagicus.txt', 'w') 
for date in date2valid_pair:
    valid_pair = date2valid_pair[date]
    line = date + '\t' + '\t'.join(valid_pair) + '\n'
    f.write(line)
f.close()    
    
           
            