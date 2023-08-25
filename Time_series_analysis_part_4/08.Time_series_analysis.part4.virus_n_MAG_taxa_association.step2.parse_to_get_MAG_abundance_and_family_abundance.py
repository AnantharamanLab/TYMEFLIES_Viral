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
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)
    
# Aim: Parse to get MAG abundance and family abundance
# Note: This script should be run under Mapping conda env: "conda activate /storage1/data11/yml_environments/ViWrap-Mapping"
# Processing the MAG presence with both coverage (>= 0.33) and breadth (>= 50%), 
# and also filtering the scaffold abundance by removing the highest 5% and lowest 5% fractions


def calculate_average_after_filtering(scf_abun_list_input):
    # Sort the scaffold abundance list in ascending order
    sorted_list = sorted(scf_abun_list_input)

    # Calculate the number of elements representing the highest 5% and lowest 5%
    n = len(sorted_list)
    top_5_percent = int(n * 0.05)
    bottom_5_percent = int(n * 0.95)

    # Filter out the highest 5% and lowest 5% fractions
    filtered_list = sorted_list[top_5_percent:bottom_5_percent]

    # Calculate the average value of the filtered list
    average_value = sum(filtered_list) / len(filtered_list)

    return average_value
    
def calculate_sum(scf_abun_list_input):
    # Calculate the sum of all scaffold abundances
    total_sum = sum(scf_abun_list_input)

    return total_sum    


# Step 1 Run CoverM for all result bam files 
coverm_cmds = []
bam_addrs = glob(f"rep_MAG_bams/*.bam")
for bam_addr in bam_addrs:
    stem_name = Path(bam_addr).stem
    parent_path = Path(bam_addr).parent    
    coverm_cmd = f"coverm contig -b {bam_addr} --min-read-percent-identity 93 -m metabat -o {parent_path}/{stem_name}.coverm_result.txt -q -t 1"
    coverm_cmds.append(coverm_cmd)
   
n = 30 # The number of parallel processes
for j in range(max(int(len(coverm_cmds)/n + 1), 1)):
    procs = [subprocess.Popen(i, shell=True, stdout=DEVNULL) for i in coverm_cmds[j*n: min((j+1)*n, len(coverm_cmds))] ]
    for p in procs:
        p.wait() 

# Step 2 Store CoverM result and parse to get MAG abundance       
## Step 2.1 Get the rep_MAG_list
rep_MAG_list = []
with open('/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.txt','r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('tymeflies'):
            tmp = line.split('\t')
            MAG, winner = tmp[5], tmp[14]
            if winner == 'TRUE':
                rep_MAG_list.append(MAG)
           
## Step 2.2 Store rep_MAG2scfs and rep_MAG2tax dict, and get rep_MAG2family and family2rep_MAGs dict
rep_MAG2scfs = {}
rep_MAG2tax = {}
with open('/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.CheckM_passed.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('IMG'):
            tmp = line.split('\t')
            MAG, scfs, tax = tmp[0], tmp[4], tmp[1]
            if MAG in rep_MAG_list:
                rep_MAG2scfs[MAG] = scfs.split(',')
                rep_MAG2tax[MAG] = tax
                
rep_MAG2family = {}
family2rep_MAGs = defaultdict(list)
for MAG in rep_MAG2tax:
    tax = rep_MAG2tax[MAG]
    family = tax.split(';')[3] + ';' + tax.split(';')[4]
    rep_MAG2family[MAG] = family
    family2rep_MAGs[family].append(MAG)
                
## Step 2.3 Get IMG2read_no dict 
IMG2read_no = {} # IMG => read_no
with open('/storage1/data11/TYMEFLIES_phage/TYMEFLIES_metagenome_info.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('IMG'):
            tmp = line.split('\t')
            IMG, read_no = tmp[0], tmp[13]
            IMG2read_no[IMG] = read_no
lines.close()           
            
## Step 2.4 Parse to get the MAG abundance
MAG2IMG2abun = defaultdict(dict) # MAG => IMG => abun (normalized)
for bam_addr in bam_addrs:
    stem_name = Path(bam_addr).stem
    parent_path = Path(bam_addr).parent   
    IMG = stem_name.split('.')[0]
    # Store scf2abun and scf2length dict
    scf2abun = {}
    scf2length = {}
    with open(f"{parent_path}/{stem_name}.coverm_result.txt", 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('contig'):
                tmp = line.split('\t')
                scf, abun, length = tmp[0], tmp[2], tmp[1]
                scf2abun[scf] = float(abun)
                scf2length[scf] = int(length)
    
    for MAG in rep_MAG2scfs:
        abun = 0 # MAG abundance
        abun_normalized = 0 # MAG abundance (normalized)
        
        total_length = 0 # The total length of scaffolds in the MAG
        total_length_not_zero = 0 # The total length of scaffold with coverage higher than 0 in the MAG
        
        scfs = rep_MAG2scfs[MAG]
        scf_abun_list = [] # Store all the scf abun
        for scf in scfs:
            scf_abun_list.append(scf2abun[scf])
            total_length += scf2length[scf]
            if scf2abun[scf] > 0:
                total_length_not_zero += scf2length[scf]
        
        if calculate_sum(scf_abun_list) > 0:
            abun = calculate_average_after_filtering(scf_abun_list) # Store the abun of the MAG after filtering the top and bottom 10% fractions
            
        if abun > 0 and float(total_length_not_zero / total_length) > 0.5: # If the non-0-coverage scaffold fraction is high than 70%
            abun_normalized = float(abun) / (float(IMG2read_no[IMG]) / 100000000) # Normalized by 100M reads
        
        if abun_normalized > 0.33:
            MAG2IMG2abun[MAG][IMG] = float(abun_normalized)
        else:
            MAG2IMG2abun[MAG][IMG] = 0
            
## Step 2.5 Write down the MAG2IMG2abun dict
#os.mkdir('virus_n_MAG_tax_association')
with open('virus_n_MAG_tax_association/MAG2IMG2abun.txt', 'w') as file:
    ## Write table header
    header = '\t'.join(sorted(IMG2read_no.keys()))
    file.write('head\t' + header + '\n')

    ## Write table rows
    for MAG in sorted(MAG2IMG2abun.keys()):
        row_values = [str(MAG2IMG2abun.get(MAG, {}).get(IMG, '')) for IMG in IMG2read_no.keys()]
        row = MAG + '\t' + '\t'.join(row_values) + '\n'
        file.write(row)        
                
# Step 3 Parse to get the family abundance
## Step 3.1 Parse to get the family abundance
family2IMG2abun = defaultdict(dict) # family => IMG => abun (normalized)
for family in family2rep_MAGs:
    for IMG in IMG2read_no:
        abun_for_the_family = 0 # The abun for the family (normalized)
        rep_MAGs = family2rep_MAGs[family]
        for MAG in rep_MAGs:        
            abun_normalized_for_the_MAG = MAG2IMG2abun[MAG][IMG]
            abun_for_the_family += abun_normalized_for_the_MAG
        family2IMG2abun[family][IMG] = abun_for_the_family           

## Step 3.2 Write down the family2IMG2abun dict
with open('virus_n_MAG_tax_association/family2IMG2abun.txt', 'w') as file:
    ## Write table header
    header = '\t'.join(sorted(IMG2read_no.keys()))
    file.write('head\t' + header + '\n')

    ## Write table rows
    for family in sorted(family2IMG2abun.keys()):
        row_values = [str(family2IMG2abun.get(family, {}).get(IMG, '')) for IMG in IMG2read_no.keys()]
        row = family + '\t' + '\t'.join(row_values) + '\n'
        file.write(row) 

