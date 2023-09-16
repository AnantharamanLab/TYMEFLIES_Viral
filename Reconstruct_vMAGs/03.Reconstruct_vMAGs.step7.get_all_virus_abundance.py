#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    warnings.filterwarnings("ignore")
    from pathlib import Path
    from glob import glob
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)
    
# Aim: Get all virus normalized abundance
     

def store_seq(input_seq_file): # The input sequence file should be a file with full path
    head = "" # Store the header line
    seq_dict = {} # Store the sequence dict
    
    with open(input_seq_file, "r") as seq_lines:
        for line in seq_lines:
            line = line.rstrip("\n") # Remove "\n" in the end
            if ">" in line:
                if (" " or "\t") in line: # Break at the first " " or "\t"
                    spliter = ""
                    for i in range(len(line)):
                        if line[i] == " " or line[i] == "\t":
                            spliter = line[i]
                            break 
                           
                    head = line.split(f'{spliter}', 1)[0]
                    seq_dict[head] = ""
                else:
                    head = line
                    seq_dict[head] = ""
            else:
                seq_dict[head] += line
            
    seq_lines.close()
    
    return seq_dict
    
def get_cov_for_scaffold_fragment(scf_fragment, IMG):                                                                            
    cov = 0 # The coverage of scf_fragment
   
    # Get all the prophage coverage results within this sample (IMG)
    propagate_result_file = f"/storage1/data11/TYMEFLIES_phage/{IMG}/PropagAtE_result/PropagAtE_result.tsv"
    prophage2cov = {} # prophage => cov
    with open(propagate_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('prophage'):
                tmp = line.split('\t')
                prophage, cov = tmp[0], tmp[7]
                prophage2cov[prophage] = cov
    lines.close()            

    if scf_fragment in prophage2cov:
        cov = prophage2cov[scf_fragment]             
    else:
        exit("The input scf_fragment is not in the dict prophage2cov")
    return cov

# Step 1 Get IMG2read_no dict
IMG2read_no = {} # IMG => read_no
with open('/storage1/data11/TYMEFLIES_phage/TYMEFLIES_metagenome_info.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        IMG, read_no = tmp[0], tmp[13]
        IMG2read_no[IMG] = read_no
lines.close()


# Step 2 Store the normalized depth for each scaffold
scf2depth_normalized = {} # scf => depth_normalized
depth_file_addrs = glob('/storage1/data11/TYMEFLIES_phage/*/*.id97.coverm_depth.txt')
for depth_file_addr in depth_file_addrs:
    IMG = Path(depth_file_addr).stem.split('.', 1)[0]
    read_no = IMG2read_no[IMG]
    with open(depth_file_addr, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('contigName'):
                tmp = line.split('\t')
                scf, depth = tmp[0], tmp[3]
                depth_normalized = float(depth) / (float(read_no) / 100000000) # Normalized by 100M reads
                scf2depth_normalized[scf] = depth_normalized
    lines.close()
    
    
# Step 3 Get all virus normalized abundance       
# Step 3.1 Store all virus to long scf dict
all_virus_fasta_addrs = []
with open('all_virus_fasta_addr.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        all_virus_fasta_addrs.append(line)
lines.close()   

viral_gn2scfs_long = {} # viral_gn => [scfs_long]
for each_virus_fasta_addr in all_virus_fasta_addrs:
    each_virus_fasta_seq = store_seq(each_virus_fasta_addr)
    viral_gn = Path(each_virus_fasta_addr).stem 
    scfs_long = []
    for header in each_virus_fasta_seq:
        scf_long = header.replace('>', '', 1)
        scfs_long.append(scf_long)
    viral_gn2scfs_long[viral_gn] = scfs_long
         
# Step 3.2 Get viral_gn2depth_normalized dict and write it down
viral_gn2depth_normalized = {} # viral_gn => depth_normalized
for viral_gn in viral_gn2scfs_long:
    scfs_long = viral_gn2scfs_long[viral_gn]
    depth_normalized = 0
    depth_normalized_sum_for_all_scfs = 0
    for scf_long in scfs_long:
        cov = 0 # Store the normalized coverage for this scf_long
        scf = scf_long.rsplit('__', 1)[1]
        IMG = scf_long.split('__', 1)[0]
        if '_fragment_' not in scf:
            cov = scf2depth_normalized[scf]
        elif '_fragment_' in scf:    
            cov = get_cov_for_scaffold_fragment(scf, IMG)
            if cov != 'NA':
                cov = float(cov) / (float(IMG2read_no[IMG]) / 100000000) # Normalized by 100M reads
            else:
                cov = 0
        depth_normalized_sum_for_all_scfs = depth_normalized_sum_for_all_scfs + cov 
    depth_normalized = depth_normalized_sum_for_all_scfs / len(scfs_long) 
    viral_gn2depth_normalized[viral_gn] = depth_normalized

f = open('viral_gn2depth_normalized.txt', 'w')
for viral_gn in viral_gn2depth_normalized:
    line = viral_gn + '\t' + str(viral_gn2depth_normalized[viral_gn]) + '\n'
    f.write(line)
f.close()           