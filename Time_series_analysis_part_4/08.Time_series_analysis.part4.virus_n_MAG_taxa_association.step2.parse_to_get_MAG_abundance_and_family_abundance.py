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
# Processing the MAG presence with breadth (>= 10%), 
  

# Step 1 Run CoverM for all result bam files 
coverm_cmds = []
bam_addrs = glob(f"rep_MAG_bams/*.bam")
for bam_addr in bam_addrs:
    stem_name = Path(bam_addr).stem
    parent_path = Path(bam_addr).parent    
    coverm_cmd = f"coverm contig -b {bam_addr} --min-read-percent-identity 93 -m metabat -o {parent_path}/{stem_name}.coverm_result.txt -q -t 1"
    coverm_cmds.append(coverm_cmd)
'''   
n = 30 # The number of parallel processes
for j in range(max(int(len(coverm_cmds)/n + 1), 1)):
    procs = [subprocess.Popen(i, shell=True, stdout=DEVNULL) for i in coverm_cmds[j*n: min((j+1)*n, len(coverm_cmds))] ]
    for p in procs:
        p.wait() 
'''
# Step 2 Store CoverM result and parse to get MAG abundance             
## Step 2.1 Store rep_MAG2scfs and rep_MAG2tax dict, and get rep_MAG2family and family2rep_MAGs dict
rep_MAG2scfs = {}
rep_MAG2tax = {}
with open('/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.rep_MAG.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('IMG'):
            tmp = line.split('\t')
            MAG, scfs, tax = tmp[0], tmp[4], tmp[1]
            rep_MAG2scfs[MAG] = scfs.split(',')
            rep_MAG2tax[MAG] = tax
                
rep_MAG2family = {}
family2rep_MAGs = defaultdict(list)
for MAG in rep_MAG2tax:
    tax = rep_MAG2tax[MAG]
    family = tax.split(';')[3] + ';' + tax.split(';')[4]
    rep_MAG2family[MAG] = family
    family2rep_MAGs[family].append(MAG)
                
## Step 2.2 Get IMG2read_no dict 
IMG2read_no = {} # IMG => read_no
with open('/storage1/data11/TYMEFLIES_phage/TYMEFLIES_metagenome_info.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('IMG'):
            tmp = line.split('\t')
            IMG, read_no = tmp[0], tmp[13]
            IMG2read_no[IMG] = read_no
lines.close()           
            
## Step 2.3 Parse to get the MAG abundance
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
        for scf in scfs:
            total_length += scf2length[scf]
            if scf2abun[scf] > 0:
                total_length_not_zero += scf2length[scf]
                   
        for scf in scfs:
            scf_length_fraction = float(scf2length[scf] / total_length) 
            abun += scf_length_fraction * scf2abun[scf]
            
        if abun > 0 and float(total_length_not_zero / total_length) > 0.1: # If the non-0-coverage scaffold fraction is high than 10%
            abun_normalized = float(abun) / (float(IMG2read_no[IMG]) / 100000000) # Normalized by 100M reads
        
        if abun_normalized > 0:
            MAG2IMG2abun[MAG][IMG] = float(abun_normalized)
        else:
            MAG2IMG2abun[MAG][IMG] = 0
            
## Step 2.4 Write down the MAG2IMG2abun dict
#os.mkdir('virus_n_MAG_tax_association')
with open('virus_n_MAG_tax_association/MAG2IMG2abun.txt', 'w') as file:
    ## Write table header
    header = '\t'.join(sorted(IMG2read_no.keys()))
    file.write('head\t' + header + '\n')

    ## Write table rows
    for MAG in sorted(MAG2IMG2abun.keys()):
        row_values = [str(MAG2IMG2abun.get(MAG, {}).get(IMG, '')) for IMG in sorted(IMG2read_no.keys())]
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
        row_values = [str(family2IMG2abun.get(family, {}).get(IMG, '')) for IMG in sorted(IMG2read_no.keys())]
        row = family + '\t' + '\t'.join(row_values) + '\n'
        file.write(row) 
        
# Step 4 Store and write down family2season2abun dict
## Step 4.1 Store season/year_season to img id info
season2num = {}  # season => num_metagenome; Store how many samples (metagenomes) are there in each season for the whole datasets
season2img_id = {}  # season => collection of img_id separated by "\t"
year_season2num = {}  # year_season => num_metagenome
year_season2img_id = {}  # year_season (for example, "2010-Ice-on") => collection of img_id separated by "\t"
year2num = {}  # year => num_metagenome
year2img_id = {}  # year (for example, "2010") => collection of img_id separated by "\t"

with open("/storage1/data11/TYMEFLIES_phage/TYMEFLIES_metagenome_info.txt") as IN:
    for line in IN:
        line = line.strip()
        if not line.startswith('IMG'):
            tmp = line.split('\t')
            img_id = tmp[0]
            date = tmp[8]
            year = date.split('-')[0]
            season = tmp[10]
            # Increment the count of metagenomes for the current season
            if season in season2num:
                season2num[season] += 1
            else:
                season2num[season] = 1

            # Append the current IMG ID to the list of IMG IDs for the current season
            if season in season2img_id:
                season2img_id[season] += "\t" + img_id
            else:
                season2img_id[season] = img_id

            # Create a year-season combination
            year_season = year + '-' + season

            # Increment the count of metagenomes for the current year-season
            if year_season in year_season2num:
                year_season2num[year_season] += 1
            else:
                year_season2num[year_season] = 1

            # Append the current IMG ID to the list of IMG IDs for the current year-season
            if year_season in year_season2img_id:
                year_season2img_id[year_season] += "\t" + img_id
            else:
                year_season2img_id[year_season] = img_id

            # Increment the count of metagenomes for the current year
            if year in year2num:
                year2num[year] += 1
            else:
                year2num[year] = 1

            # Append the current IMG ID to the list of IMG IDs for the current year
            if year in year2img_id:
                year2img_id[year] += "\t" + img_id
            else:
                year2img_id[year] = img_id

# Step 4.2 Get family abundances for each season (the number of metagenome each season was normalized)
family2season2abun = defaultdict(dict)  # family => season => abun (normalized abundance)
for family in sorted(family2IMG2abun):
    for season in sorted(season2num):
        num_metagenome = season2num[season]  # num of metagenomes in this season
        abun_sum_for_season = 0  # Store the sum of family abun from all metagenomes from this season
        abun_for_season_normalized = 0  # Store the normalized family abun for this season

        IMG_ID = season2img_id[season].split('\t')
        for img_id in IMG_ID:
            abun = family2IMG2abun[family].get(img_id, 0)
            if abun:
                abun_sum_for_season += abun

        if abun_sum_for_season:
            abun_for_season_normalized = abun_sum_for_season / num_metagenome

        family2season2abun[family][season] = abun_for_season_normalized

# Step 4.3 Write down Family abundances (normalized by read number per metagenome and metagenome number per season) for each season
Season = ["Spring", "Clearwater", "Early Summer", "Late Summer", "Fall", "Ice-on"]
with open("virus_n_MAG_tax_association/family2season2abun.txt", 'w') as OUT:
    row = "\t".join(Season)
    OUT.write("Head\t" + row + "\n")
    for tmp1 in sorted(family2season2abun):
        OUT.write(tmp1 + "\t")
        tmp = []
        for tmp2 in Season:
            if tmp2 in family2season2abun[tmp1]:
                tmp.append(str(family2season2abun[tmp1][tmp2]))
            else:
                tmp.append("0")
        OUT.write("\t".join(tmp) + "\n")


