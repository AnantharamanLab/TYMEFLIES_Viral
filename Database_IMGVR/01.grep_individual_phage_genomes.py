#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    warnings.filterwarnings("ignore")
    from pathlib import Path
    from collections import defaultdict
    from glob import glob    
    
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)
    
#########################################################
# Step 1. Check phage genome name congruency            #
# Step 2. Write down "IMGVR_all_genome_scaffolds.txt"   #
#########################################################

# Functions
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

# Step 1. Check phage genome name congruency      
## Step 1.1 Determine the genome name from IMGVR_all_nucleotides_header file    
Gn2scfs = defaultdict(list) # gn => [scfs]
with open("IMGVR_all_nucleotides_header.txt", "r") as f:
    for line in f:
        scf = line.strip()
        gn = ""
        
        # Situation 1. Scf contains UViG
        # All scf has 2 or 3 "|", the content before the 1st "|" is the genome name
        if "UViG" in scf:
            gn = scf.split("|")[0]
        else:
            # Sitation 2. Scf has IMG ID inside, and before IMG ID, the content is not IMG ID again, the content is the genome name
            if scf.startswith(".+?|\d\d\d\d\d\d\d\d\d\d|") and not scf.startswith("\d\d\d\d\d\d\d\d\d\d|\d\d\d\d\d\d\d\d\d\d|"):
                gn = scf.split("|")[0]
            elif scf.startswith("\d\d\d\d\d\d\d\d\d\d|"): # Sitation 3. Scf has IMG ID inside and in the front, and it is a single-scaffold phage name, the phage genome name is the same as the scaffold name
                gn = scf
            elif "|" not in scf: # Situation 4. Scf doesn't have "|" inside, the whole scf name is the phage genome name
                gn = scf
            else: # Other situations, manually check
                gn = scf

        if gn:
            Gn2scfs[gn].append(scf)
f.close()  

## Step 1.2 Store the genome name from "IMGVR_all_Sequence_information.tsv", and check if our genome name is the same as this
Gn_from_all_seq_info = {} # gn => 1
with open("IMGVR_all_Sequence_information.tsv", "r") as f:
    for line in f:
        if not line.startswith("UVIG"):
            gn = line.strip().split("\t")[0]
            Gn_from_all_seq_info[gn] = 1

# Check if our gn ids are all present in theirs
for gn1 in sorted(Gn2scfs.keys()):
    if gn1 not in Gn_from_all_seq_info:
        print(f"{gn1} is not in IMGVR_all_Sequence_information.tsv")

# Check if their gn ids are all present in ours
for gn2 in sorted(Gn_from_all_seq_info.keys()):
    if gn2 not in Gn2scfs:
        print(f"{gn2} is not in our Gn2scfs dict")

# Check the numbers of elements in both dict
element_num_Gn2scfs = len(Gn2scfs)
element_num_Gn_from_all_seq = len(Gn_from_all_seq_info)     
print(f"Our dict contains {element_num_Gn2scfs} elements\n")
print(f"Their dict contains {element_num_Gn_from_all_seq} elements\n")  


# Step 2. Write down "IMGVR_all_genome_scaffolds.txt" to store scaffolds to genome map
scaffold_no = 0 # Summarize all the scaffold nums in all genomes
with open("IMGVR_all_genome_scaffolds.txt", "w") as f:
    f.write("viral genome\tscaffold number\tscaffolds\n")
    for gn in sorted(Gn2scfs.keys()):
        scf_num = len(Gn2scfs[gn])
        scaffold_no = scaffold_no + scf_num
        f.write(f"{gn}\t{scf_num}\t{','.join(Gn2scfs[gn])}\n")
f.close() 
       
print(f"All the scaffold num is {scaffold_no}\n")       