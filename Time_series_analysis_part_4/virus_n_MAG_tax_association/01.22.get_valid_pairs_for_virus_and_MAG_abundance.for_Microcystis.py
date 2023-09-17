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
# For Microcystis:
# (1) Microcystis_MAG (2) psbA_containing_Microcystis_virus (3) no_psbA_containing_Microcystis_virus


# Step 1 Store Microcystis_MAG abundance from interpolation results
## Step 1.1 Store Microcystis_MAG_year2highest_y_value dict
Microcystis_MAG_year2highest_y_value = {} # year => highest_y_value
with open('Microcystis_virus_n_MAG/Microcystis_MAG.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        Microcystis_MAG_year2highest_y_value[year] = highest_y_value

## Step 1.2 Store Microcystis_MAG_date2y_value dict
Microcystis_MAG_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
Microcystis_MAG_interpolation_result_files = glob('Microcystis_virus_n_MAG/*.Microcystis_MAG.interpolation_results.txt')
for Microcystis_MAG_interpolation_result_file in Microcystis_MAG_interpolation_result_files:
    year = Microcystis_MAG_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(Microcystis_MAG_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            Microcystis_MAG_date2y_value[f"{year}|{date}"] = y_value
            
  
# Step 2 Store psbA_containing_Microcystis_viral_gn abundance from interpolation results
## Step 2.1 Store psbA_containing_Microcystis_viral_gn_year2highest_y_value dict
psbA_containing_Microcystis_viral_gn_year2highest_y_value = {} # year => highest_y_value
with open('Microcystis_virus_n_MAG/psbA_containing_Microcystis_viral_gn.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        psbA_containing_Microcystis_viral_gn_year2highest_y_value[year] = highest_y_value

## Step 2.2 Store psbA_containing_Microcystis_viral_gn_date2y_value dict
psbA_containing_Microcystis_viral_gn_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
psbA_containing_Microcystis_viral_gn_interpolation_result_files = glob('Microcystis_virus_n_MAG/*.psbA_containing_Microcystis_viral_gn.interpolation_results.txt')
for psbA_containing_Microcystis_viral_gn_interpolation_result_file in psbA_containing_Microcystis_viral_gn_interpolation_result_files:
    year = psbA_containing_Microcystis_viral_gn_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(psbA_containing_Microcystis_viral_gn_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            psbA_containing_Microcystis_viral_gn_date2y_value[f"{year}|{date}"] = y_value
            

# Step 3 Store no_psbA_containing_Microcystis_viral_gn abundance from interpolation results
## Step 3.1 Store no_psbA_containing_Microcystis_viral_gn_year2highest_y_value dict
no_psbA_containing_Microcystis_viral_gn_year2highest_y_value = {} # year => highest_y_value
with open('Microcystis_virus_n_MAG/no_psbA_containing_Microcystis_viral_gn.year2highest_y_value.txt') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        year, highest_y_value = tmp[0], tmp[1]
        no_psbA_containing_Microcystis_viral_gn_year2highest_y_value[year] = highest_y_value

## Step 3.2 Store no_psbA_containing_Microcystis_viral_gn_date2y_value dict
no_psbA_containing_Microcystis_viral_gn_date2y_value = {} # year_date (for instance, '2010|-45') => y_value
no_psbA_containing_Microcystis_viral_gn_interpolation_result_files = glob('Microcystis_virus_n_MAG/*.no_psbA_containing_Microcystis_viral_gn.interpolation_results.txt')
for no_psbA_containing_Microcystis_viral_gn_interpolation_result_file in no_psbA_containing_Microcystis_viral_gn_interpolation_result_files:
    year = no_psbA_containing_Microcystis_viral_gn_interpolation_result_file.split('/')[1].split('.')[0] 
    with open(no_psbA_containing_Microcystis_viral_gn_interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            date, y_value = tmp[0], tmp[1]
            no_psbA_containing_Microcystis_viral_gn_date2y_value[f"{year}|{date}"] = y_value

            
# Step 4 Filter valid pairs
date2valid_pair = {} # date => [Microcystis_MAG_abun, psbA_containing_Microcystis_viral_gn_abun, no_psbA_containing_Microcystis_viral_gn_abun]
for date in Microcystis_MAG_date2y_value:
    Microcystis_MAG_abun_ratio = Microcystis_MAG_date2y_value[date]
    psbA_containing_Microcystis_viral_gn_abun_ratio = psbA_containing_Microcystis_viral_gn_date2y_value[date]
    no_psbA_containing_Microcystis_viral_gn_abun_ratio = no_psbA_containing_Microcystis_viral_gn_date2y_value[date]
    
    year = date.split('|')[0]    
    if Microcystis_MAG_abun_ratio != 'nan' and psbA_containing_Microcystis_viral_gn_abun_ratio != 'nan' and no_psbA_containing_Microcystis_viral_gn_abun_ratio != 'nan':
        if Microcystis_MAG_year2highest_y_value[year] != 'nan' and psbA_containing_Microcystis_viral_gn_year2highest_y_value[year] != 'nan' and no_psbA_containing_Microcystis_viral_gn_year2highest_y_value[year] != 'nan': 
            Microcystis_MAG_abun = float(Microcystis_MAG_abun_ratio) * float(Microcystis_MAG_year2highest_y_value[year]) / 100    
            psbA_containing_Microcystis_viral_gn_abun = float(psbA_containing_Microcystis_viral_gn_abun_ratio) * float(psbA_containing_Microcystis_viral_gn_year2highest_y_value[year]) / 100 
            no_psbA_containing_Microcystis_viral_gn_abun = float(no_psbA_containing_Microcystis_viral_gn_abun_ratio) * float(no_psbA_containing_Microcystis_viral_gn_year2highest_y_value[year]) / 100                     
            if float(Microcystis_MAG_abun_ratio) > 10 and float(psbA_containing_Microcystis_viral_gn_abun_ratio) > 10 and float(no_psbA_containing_Microcystis_viral_gn_abun_ratio) > 10:
                date2valid_pair[date] = [str(Microcystis_MAG_abun), str(psbA_containing_Microcystis_viral_gn_abun), str(no_psbA_containing_Microcystis_viral_gn_abun)]
            
f = open('Microcystis_virus_n_MAG/psbA_date2valid_pair.for_Microcystis.txt', 'w') 
for date in date2valid_pair:
    valid_pair = date2valid_pair[date]
    line = date + '\t' + '\t'.join(valid_pair) + '\n'
    f.write(line)
f.close()    
    
           
            