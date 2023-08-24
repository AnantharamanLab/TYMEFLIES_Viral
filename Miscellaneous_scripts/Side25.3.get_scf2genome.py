#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    import statistics
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)

# Aim: Get Scf2genome.stb from "All_phage_species_rep_gn.fasta"


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
    
    
# Step 1 Store the fasta seq
All_phage_species_rep_gn_seq = store_seq('reference_fasta_for_metapop/All_phage_species_rep_gn.fasta')


# Step 2 Parse to get scf2genome dict
scf2genome = {} # scf => genome
for header in All_phage_species_rep_gn_seq:
    scf = header.replace('>', '', 1)
    genome = scf.rsplit('__', 1)[0]
    scf2genome[scf] = genome
    

# Step 3 Write down Scf2genome.stb
f = open('Summer_vs_Winter_Fst_analysis/Scf2genome.stb', 'w')
for scf in scf2genome:
    genome = scf2genome[scf]
    line = scf + '\t' + genome
    f.write(line + '\n')
f.close()    


