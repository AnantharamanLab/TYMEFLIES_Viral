#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    from collections import defaultdict
    from pathlib import Path
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: Replace the bin seq headers of Robin's MAGs


def store_seq_with_full_head(input_seq_file): # The input sequence file should be a file with full path
    head = "" # Store the header line
    seq_dict = {} # Store the sequence dict
    
    with open(input_seq_file, "r") as seq_lines:
        for line in seq_lines:
            line = line.replace("\n","") # Remove "\n"
            if line[0] == ">":
                head = line
                seq_dict[head] = ""
            else:
                seq_dict[head] += line
            
    seq_lines.close()
    
    return seq_dict
    
def write_down_seq(seq_dict, path_to_file): 
    # Two inputs are required:
    # (1) The dict of the sequence
    # (2) The path that you want to write your sequence down
    
    seq_file = open(path_to_file,"w")
    for head in seq_dict:
        seq_file.write(head + "\n")
        seq_file.write(seq_dict[head] + "\n")
    seq_file.close()      
    
# Step 1 Store the Robin_scf_name to IMG_scf_name map
Robin_scf_name2IMG_scf_name = {} # An example: ME2005-08-02_3300034104_group2_bin150_scaffold_7774_c1 => Ga0453673_0007353
with open('bin_info/bin_scaffold_key.csv', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('Assembly'):
            tmp = line.split(',')
            Robin_scf_name, IMG_scf_name = tmp[2], tmp[1]
            Robin_scf_name2IMG_scf_name[Robin_scf_name] = IMG_scf_name
lines.close()


# Step 2 Store each fasta file and replace the headers
All_fna_addrs = glob('ME*/*.fna')
for each_fna_addr in All_fna_addrs:
    each_fna_seq = store_seq_with_full_head(each_fna_addr)
    each_fna_seq_new = {}
    for header in each_fna_seq:
        header_wo_arrow = header.replace('>', '', 1)
        if header_wo_arrow in Robin_scf_name2IMG_scf_name:
            header_new = '>' + Robin_scf_name2IMG_scf_name[header_wo_arrow]
            each_fna_seq_new[header_new] = each_fna_seq[header]
        else:
            print(f"{header_wo_arrow} is not in bin_scaffold_key.csv")
    IMG = Path(each_fna_addr).stem.split('_')[1]
    each_fna_addr_new = Path(IMG) / str(Path(each_fna_addr).stem + '.fasta')
    if not os.path.exists(IMG):
        os.mkdir(IMG)
    write_down_seq(each_fna_seq_new, each_fna_addr_new)               