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
# For Planktophila:
# (1) Planktophila_MAG (2) katG_containing_Planktophila_viral_gn (3) no_katG_containing_Planktophila_viral_gn


# Step 1 Store Planktophila_MAG abundance from interpolation results
## Step 1.1 Store Planktophila_MAG_year2highest_y_value dict
Planktophila_MAG_year2highest_y_value = {} # year => highest_y_value
with open('Planktophila_virus_n_MAG/Planktophila_MAG.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        Planktophila_MAG_year2highest_y_value[year] = highest_y_value

## Step 1.2 Store Planktophila_MAG_date2y_value dict
Planktophila_MAG_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
Planktophila_MAG_interpolation_result_files = glob('Planktophila_virus_n_MAG/*.Planktophila_MAG.interpolation_results.txt')
for Planktophila_MAG_interpolation_result_file in Planktophila_MAG_interpolation_result_files:
    year = Planktophila_MAG_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(Planktophila_MAG_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            Planktophila_MAG_date2y_value[f"{year}|{date}"] = y_value
            
  
# Step 2 Store katG_containing_Planktophila_viral_gn abundance from interpolation results
## Step 2.1 Store katG_containing_Planktophila_viral_gn_year2highest_y_value dict
katG_containing_Planktophila_viral_gn_year2highest_y_value = {} # year => highest_y_value
with open('Planktophila_virus_n_MAG/katG_containing_Planktophila_viral_gn.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        katG_containing_Planktophila_viral_gn_year2highest_y_value[year] = highest_y_value

## Step 2.2 Store katG_containing_Planktophila_viral_gn_date2y_value dict
katG_containing_Planktophila_viral_gn_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
katG_containing_Planktophila_viral_gn_interpolation_result_files = glob('Planktophila_virus_n_MAG/*.katG_containing_Planktophila_viral_gn.interpolation_results.txt')
for katG_containing_Planktophila_viral_gn_interpolation_result_file in katG_containing_Planktophila_viral_gn_interpolation_result_files:
    year = katG_containing_Planktophila_viral_gn_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(katG_containing_Planktophila_viral_gn_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            katG_containing_Planktophila_viral_gn_date2y_value[f"{year}|{date}"] = y_value
            

# Step 3 Store no_katG_containing_Planktophila_viral_gn abundance from interpolation results
## Step 3.1 Store no_katG_containing_Planktophila_viral_gn_year2highest_y_value dict
no_katG_containing_Planktophila_viral_gn_year2highest_y_value = {} # year => highest_y_value
with open('Planktophila_virus_n_MAG/no_katG_containing_Planktophila_viral_gn.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        no_katG_containing_Planktophila_viral_gn_year2highest_y_value[year] = highest_y_value

## Step 3.2 Store no_katG_containing_Planktophila_viral_gn_date2y_value dict
no_katG_containing_Planktophila_viral_gn_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
no_katG_containing_Planktophila_viral_gn_interpolation_result_files = glob('Planktophila_virus_n_MAG/*.no_katG_containing_Planktophila_viral_gn.interpolation_results.txt')
for no_katG_containing_Planktophila_viral_gn_interpolation_result_file in no_katG_containing_Planktophila_viral_gn_interpolation_result_files:
    year = no_katG_containing_Planktophila_viral_gn_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(no_katG_containing_Planktophila_viral_gn_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            no_katG_containing_Planktophila_viral_gn_date2y_value[f"{year}|{date}"] = y_value

            
# Step 4 Filter valid pairs
date2valid_pair = {} # date => [Planktophila_MAG_abun, katG_containing_Planktophila_viral_gn_abun, no_katG_containing_Planktophila_viral_gn_abun]
for date in Planktophila_MAG_date2y_value:
    Planktophila_MAG_abun_ratio = Planktophila_MAG_date2y_value[date]
    katG_containing_Planktophila_viral_gn_abun_ratio = katG_containing_Planktophila_viral_gn_date2y_value[date]
    no_katG_containing_Planktophila_viral_gn_abun_ratio = no_katG_containing_Planktophila_viral_gn_date2y_value[date]
    
    year = date.split('|')[0]    
    if Planktophila_MAG_abun_ratio != 'nan' and katG_containing_Planktophila_viral_gn_abun_ratio != 'nan' and no_katG_containing_Planktophila_viral_gn_abun_ratio != 'nan':
        if Planktophila_MAG_year2highest_y_value[year] != 'nan' and katG_containing_Planktophila_viral_gn_year2highest_y_value[year] != 'nan' and no_katG_containing_Planktophila_viral_gn_year2highest_y_value[year] != 'nan': 
            Planktophila_MAG_abun = float(Planktophila_MAG_abun_ratio) * float(Planktophila_MAG_year2highest_y_value[year]) / 100    
            katG_containing_Planktophila_viral_gn_abun = float(katG_containing_Planktophila_viral_gn_abun_ratio) * float(katG_containing_Planktophila_viral_gn_year2highest_y_value[year]) / 100 
            no_katG_containing_Planktophila_viral_gn_abun = float(no_katG_containing_Planktophila_viral_gn_abun_ratio) * float(no_katG_containing_Planktophila_viral_gn_year2highest_y_value[year]) / 100                     
            if float(Planktophila_MAG_abun_ratio) > 10 and float(katG_containing_Planktophila_viral_gn_abun_ratio) > 10 and float(no_katG_containing_Planktophila_viral_gn_abun_ratio) > 10:
                date2valid_pair[date] = [str(Planktophila_MAG_abun), str(katG_containing_Planktophila_viral_gn_abun), str(no_katG_containing_Planktophila_viral_gn_abun)]
            
f = open('Planktophila_virus_n_MAG/katG_date2valid_pair.for_Planktophila.txt', 'w') 
for date in date2valid_pair:
    valid_pair = date2valid_pair[date]
    line = date + '\t' + '\t'.join(valid_pair) + '\n'
    f.write(line)
f.close()    
    
           
            