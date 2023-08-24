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
    
    
# Aim: Run GTDB-Tk for Robin's MAGs
# Note: This script should be run under conda env "gtdbtk_2.1.1"


# Step 1 Copy all the MAGs into a new folder 
all_fasta_addrs = glob('33*/*.fasta')
os.mkdir('all_MAGs_fasta')
for each_fasta_addr in all_fasta_addrs:
    each_fasta_addr_new =  'all_MAGs_fasta' + '/' + str(Path(each_fasta_addr).stem + '.fasta')
    os.system(f"cp {each_fasta_addr} {each_fasta_addr_new}")
    
    
# Step 2 Run GTDB-Tk
gtdbtk_cmd = 'gtdbtk classify_wf --genome_dir all_MAGs_fasta --out_dir all_MAGs_fasta_gtdbtk_output_dir -x .fasta --cpus 15'
os.system(gtdbtk_cmd)         