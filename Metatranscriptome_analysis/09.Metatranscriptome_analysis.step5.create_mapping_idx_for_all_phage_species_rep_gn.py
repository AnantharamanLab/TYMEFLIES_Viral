#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    from pathlib import Path    
    import pyfastx # For fastq and fasta reading and parsing   
    import pysam    
    warnings.filterwarnings('ignore')
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
# Aim: Making mapping idx for all_phage_species_rep_gn for 2020 metaT datasets and the 24hrs metaT datasets

def make_bowtie2_idx(fasta, working_dir, num_threads):
    file_name = Path(fasta).stem
    index_name = file_name + ".bowtie2_idx"
    
    fa = pyfastx.Fasta(fasta)
    fasta_size = fa.size
    
    if fasta_size <= 4000000000:
        # Indexing the reference sequence 
        indexing_cmd = f'bowtie2-build {fasta} {working_dir}/{index_name} --threads {num_threads} --quiet 1> /dev/null'
        os.system(indexing_cmd)
    else:
        # Indexing the reference sequence 
        indexing_cmd = f'bowtie2-build --large-index {fasta} {working_dir}/{index_name} --threads {num_threads} --quiet 1> /dev/null'
        os.system(indexing_cmd)
 
mapping_ref_all_phage_species_rep = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_24hrs_metaT_ref.based_on_all_phage_species_rep.fa'
num_threads = 10
working_dir = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome'
make_bowtie2_idx(mapping_ref_all_phage_species_rep, working_dir, num_threads)