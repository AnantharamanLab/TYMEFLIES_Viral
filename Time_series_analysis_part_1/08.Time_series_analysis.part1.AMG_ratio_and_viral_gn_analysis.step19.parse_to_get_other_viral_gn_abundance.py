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
    from statistics import mean    
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)
    
# Aim: Parse to get other viral gn abundance (viruses that do not contain any AMGs)
# Screen viral gn with the following two criteria:
# (1) Screen viral gn with its any scaffold with < 0.01 coverage 
# (2) Process the viral gn presence with both coverage (>= 0.33) and breadth (>= 50%)


def store_seq(input_seq_file): # The input sequence file should be a file with full path
    head = "" # Store the header line
    seq_dict = {} # Store the sequence dict
    
    with open(input_seq_file, "r") as seq_lines:
        for line in seq_lines:
            line = line.rstrip("\n") # Remove "\n" in the end
            if ">" in line:
                head = line.split(None, 1)[0] # Cut at the first " " or "\t", use the first part
                seq_dict[head] = ""                 
            else:
                seq_dict[head] += line
            
    seq_lines.close()
    
    return seq_dict
    
    
# Step 1 Store viral gn 2 scfs dict
## Step 1.1 Store AMG KO information
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
            
## Step 1.2 Store viral_gn2scfs dict and scf2length dict
viral_gn2scfs = defaultdict(list) # viral_gn => [scfs]
scf2length = {} # scf => length
All_phage_species_rep_gn_seq_dict = store_seq("reference_fasta_for_metapop/All_phage_species_rep_gn.fasta") 
for header in All_phage_species_rep_gn_seq_dict:
    scf = header.replace('>', '', 1)
    seq = All_phage_species_rep_gn_seq_dict[header]
    scf2length[scf] = len(seq)
    
    viral_gn = scf.split('__Ga')[0]
    if viral_gn not in AMG_containing_viral_gn: # Only keep viruses that do not contain any AMGs
        viral_gn2scfs[viral_gn].append(scf)
        
        
# Step 2 Get scf2IMG2breadth and scf2IMG2cov dicts
scf2IMG2breadth = defaultdict(dict)
scf2IMG2cov = defaultdict(dict)
breadth_and_depth_file_addrs = glob("/storage1/data11/TYMEFLIES_phage/MetaPop/03.Breadth_and_Depth/*breadth_and_depth.tsv")
for breadth_and_depth_file_addr in breadth_and_depth_file_addrs:
    IMG = Path(breadth_and_depth_file_addr).stem.split('.')[0]
    with open(breadth_and_depth_file_addr, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            scf, breadth, cov = tmp[0], tmp[2], tmp[3]
            viral_gn = scf.split('__Ga')[0]
            if viral_gn in viral_gn2scfs:
                scf2IMG2breadth[scf][IMG] = float(breadth)
                scf2IMG2cov[scf][IMG] = float(cov)


# Step 3 Get viral_gn2IMG2cov_norm_filtered (screen viral gn with the two criteria)
## Step 3.1 Store IMG2read_num dict
IMG2read_num = {}  # img -> read_num (int)
with open("Read_count_file_for_metapop.txt", "r") as infile:
    for line in infile:
        line = line.strip()
        tmp = line.split('\t')
        img_str = tmp[0].split('.')[0]
        if img_str.isdigit():
            img = img_str
            read_num = int(tmp[1])
            IMG2read_num[img] = read_num
            
## Step 3.2 Store viral_gn2IMG2breadth dict 
viral_gn2IMG2breadth = defaultdict(dict)  # viral_gn => IMG => breadth
for viral_gn in viral_gn2scfs:
    for IMG in IMG2read_num:
        scfs = viral_gn2scfs[viral_gn]
        total_breadth_multiplied_by_length = 0
        total_length = 0
        for scf in scfs:
            breadth = 0
            if scf in scf2IMG2breadth and IMG in scf2IMG2breadth[scf]:
                breadth = scf2IMG2breadth[scf][IMG]
            length = scf2length[scf]
            total_breadth_multiplied_by_length += breadth * length
            total_length += length

        breadth_for_viral_gn = 0
        if total_length:
            breadth_for_viral_gn = total_breadth_multiplied_by_length / total_length
        viral_gn2IMG2breadth[viral_gn][IMG] = breadth_for_viral_gn    
            
## Step 3.3 Store viral_gn2IMG2cov_norm_filtered dict
viral_gn2IMG2cov_norm_filtered = defaultdict(dict)
for viral_gn in viral_gn2scfs:    
    for IMG in IMG2read_num:
        scfs = viral_gn2scfs[viral_gn]
        scf_cov_norm_collection = [] # Store all the cov_norm values
        for scf in scfs:
            scf_cov = 0
            if scf in scf2IMG2cov and IMG in scf2IMG2cov[scf]:
                scf_cov = scf2IMG2cov[scf][IMG]
            scf_cov_norm = scf_cov * (100000000 / IMG2read_num[IMG])  # Normalized cov, normalized by 100M reads per metagenome
            scf_cov_norm_collection.append(scf_cov_norm)      

        # Filter viral_gn if any scaffold coverage is < 0.01
        logic = True  # Preset all scaffold coverage is >= 0.01
        for key in scf_cov_norm_collection:
            if key < 0.01:
                logic = False

        if logic:
            scf_cov_norm_mean = mean(scf_cov_norm_collection)
            breadth = 0 # the breadth of the viral gn
            if viral_gn in viral_gn2IMG2breadth and IMG in viral_gn2IMG2breadth[viral_gn]: 
                breadth = viral_gn2IMG2breadth[viral_gn][IMG]
            
            if scf_cov_norm_mean >= 0.33 and breadth >= 50:
                viral_gn2IMG2cov_norm_filtered[viral_gn][IMG] = scf_cov_norm_mean
                
                
# Step 4 Write down the result
# First, collect all unique IMGs in sorted order
all_IMGs = sorted(IMG2read_num.keys())

# Header line
header = "Head\t" + "\t".join(all_IMGs)

# Initialize the table content
table_content = [header]

# Loop through the viral_gn2IMG2cov_norm_filtered dictionary
for viral_gn in sorted(viral_gn2IMG2cov_norm_filtered.keys()):
    line = viral_gn
    for IMG in all_IMGs:
        # Check if the viral_gn and IMG exist in the dictionary
        if viral_gn in viral_gn2IMG2cov_norm_filtered and IMG in viral_gn2IMG2cov_norm_filtered[viral_gn]:
            value = viral_gn2IMG2cov_norm_filtered[viral_gn][IMG]
            # If the value is empty, replace it with '0'
            line += "\t" + str(value) if value is not None else "\t0"
        else:
            # If viral_gn or IMG does not exist, assign '0' to the cell
            line += "\t0"
    
    # Append the line to the table content
    table_content.append(line)

# Write the table to the file
with open("MetaPop/no_AMG_viral_gn2IMG2cov_norm_filtered.txt", "w") as outfile:
    outfile.write("\n".join(table_content))                
            