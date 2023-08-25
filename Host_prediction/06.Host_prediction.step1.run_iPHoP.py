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
    
# Aim: Run iPHoP to get host prediction result  
# Note: (1) This script should be run under the GTDBTk conda env for Step 1: "conda activate /storage2/scratch/zzhou388/ViWrap/conda_envs/ViWrap-GTDBTk"
#       (2) This script should be run under the iPHoP conda env for Step 2: "conda activate /storage2/scratch/zzhou388/ViWrap/conda_envs/ViWrap-iPHoP"
#       (3) This script should be run under the iPHoP conda env for Step 2: "conda activate /storage2/scratch/zzhou388/ViWrap/conda_envs/ViWrap-iPHoP"


threads = 80
custom_MAGs_dir = 'rep_MAGs'

# Step 1 Add custom MAGs to host db - make gtdbtk results
os.mkdir('custom_MAGs_GTDB-tk_results')
gtdb_cmd_1 = f'gtdbtk de_novo_wf --genome_dir {custom_MAGs_dir} --bacteria --outgroup_taxon p__Patescibacteria --out_dir custom_MAGs_GTDB-tk_results --cpus {threads} --force --extension fasta 1> /dev/null'
gtdb_cmd_2 = f'gtdbtk de_novo_wf --genome_dir {custom_MAGs_dir} --archaea --outgroup_taxon p__Altiarchaeota --out_dir custom_MAGs_GTDB-tk_results --cpus {threads} --force --extension fasta 1> /dev/null'
os.system(gtdb_cmd_1)   
os.system(gtdb_cmd_2)  


# Step 2 Add db
add_to_db_cmd = f"iphop add_to_db --fna_dir {custom_MAGs_dir} --gtdb_dir custom_MAGs_GTDB-tk_results/ --out_dir iPHoP_db_w_TYMEFLIES_rep_MAGs --db_dir /storage2/scratch/zzhou388/ViWrap/ViWrap_db/iPHoP_db/iPHoP_db/ -t {threads}"
os.system(add_to_db_cmd) 


# Step 3 Run iPHoP prediction
## Cat all Nlinked virus genomes
## all_Nlinked_virus_genomes.fasta is the concatenated virus genomes with N-linked sequences to make temporary ‘single-contig’ viral genomes
## Note: The input fasta file can be split into multiple small files and run separately to avoid memory leaks
predict_cmd = f"iphop predict --fa_file all_Nlinked_virus_genomes.fasta --db_dir iPHoP_db_w_TYMEFLIES_rep_MAGs --out_dir iphop_output/ --no_qc -t {threads}"
os.system(predict_cmd)
