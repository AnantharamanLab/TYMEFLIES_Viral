#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    import datetime
    from collections import defaultdict
    import shutil
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: Use geNomad to annotate all viral scaffolds to get taxonomy
# This script should be run under the conda env "/storage2/scratch/kosmopoulos/miniconda3/envs/genomad":
# by "conda activate /storage2/scratch/kosmopoulos/miniconda3/envs/genomad" 



def Nlinker(infolder, outdir, extension, n):
# Copied from vRhyme auxiliary scripts by Kristopher Kieft, UW-Madison
# Using n of Ns to link scaffolds
    files = os.listdir(infolder)
    len_ext = len(extension)
    files = [i for i in files if i[-len_ext:] == extension]

    if len(files) == 0:
        sys.stderr.write("\nError: No input files were identified. Verify that the input folder and extension are correct. Exiting.\n\n")
        exit()

    N_string = ''.join("N" * n)
    for f in files:
        file = infolder + "/" + f
        base = f.rsplit(".",1)[0]
        with open(file, 'r') as fasta:
            seq = ''
            for line in fasta:
                if not line.startswith(">"):
                    seq += line.strip("\n")
                    seq += N_string

        with open(outdir + "/" + base + '.linked.' + extension, "w") as outfile:
            outfile.write(">" + base + "\n" + seq[:-n] + "\n") # -n to remove last n added characters 
   
   
            
# N link fasta files 
infolder = 'all_virus_genome_fastas'
outdir = 'all_virus_genome_fastas_Nlinked'
extension = 'fasta'
n = 1000
if not os.path.exists(outdir):
    os.makedirs(outdir)
Nlinker(infolder, outdir, extension, n)

# Directory containing the FASTA files
directory = "all_virus_genome_fastas_Nlinked"

# Output file name
output_file = "all_virus_genome_fastas_Nlinked.fasta"

# Open the output file in write mode
with open(output_file, "w") as outfile:
    # Iterate over each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".fasta"):
            filepath = os.path.join(directory, filename)
            # Open each FASTA file and read its contents
            with open(filepath, "r") as infile:
                # Write the contents to the output file
                outfile.write(infile.read())
            outfile.write("\n")  # Add a newline after each file

threads = 30
genomad_db = "/storage2/scratch/kosmopoulos/databases/genomad_db"

cmd = f"genomad annotate all_virus_genome_fastas_Nlinked.fasta Taxonomic_classification/genomad_output {genomad_db} -t {threads} -q"
os.system(cmd)

