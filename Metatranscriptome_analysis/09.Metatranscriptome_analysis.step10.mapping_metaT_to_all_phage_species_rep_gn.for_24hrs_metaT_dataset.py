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
    
# Aim: Mapping metaT to all_phage_species_rep_gn for 24hrs metaT datasets (2015/8/20-21)

 
def run_bowtie2(fasta, input_read_pair, mapping_index, working_dir, sam_name, num_threads):
    # Mapping 
    mapping_cmd = f'bowtie2 -x {mapping_index} -1 {input_read_pair.split(",")[0]} -2 {input_read_pair.split(",")[1]} -S {working_dir}/{sam_name}.sam -p {int(num_threads)} --no-unal --quiet --mm 1> /dev/null'  
    os.system(mapping_cmd) 

def convert_sam_to_sorted_bam(input_sam_file, num_threads):
    # Open the SAM file in reading mode
    samfile = pysam.AlignmentFile(input_sam_file, "r")

    # Create a BAM file in writing mode
    out_bam_file = input_sam_file.replace('.sam', '.bam', 1)
    bamfile = pysam.AlignmentFile(out_bam_file, "wb", template=samfile)

    # Iterate through the records in the SAM file and write them to the BAM file
    for record in samfile:
        bamfile.write(record)

    # Close the files
    samfile.close()
    bamfile.close()

    # Sort the BAM file
    out_sorted_bam_file = input_sam_file.replace('.sam', '.sorted.bam', 1)
    pysam.sort("-o", out_sorted_bam_file, out_bam_file) 

def filter_sorted_bam(out_sorted_bam_file, filtered_bam_file, reads_mapping_identity_cutoff, aligned_length, threads):   
    reads_mapping_identity_cutoff = int(float(reads_mapping_identity_cutoff) * 100)
    threads = int(threads)
    filter_cmd = f'coverm filter --bam-files {out_sorted_bam_file} --output-bam-files {filtered_bam_file} --min-read-aligned-length {aligned_length} --min-read-percent-identity {reads_mapping_identity_cutoff} --threads {threads}'
    os.system(filter_cmd)    
    
def mapping_metaT_reads(mapping_ref, mapping_index, metaT_reads, mapping_result_dir, reads_mapping_identity_cutoff, threads):
    threads = int(threads)
    
    # Step 1 Run Bowtie2
    os.mkdir(mapping_result_dir)
    metaT_reads_list = metaT_reads.split(',')
    sam_names = []
    if len(metaT_reads_list) / 2 == 1:
        sam_name = Path(metaT_reads_list[0]).stem.rsplit('_', 1)[0]
        sam_names.append(sam_name)
        each_metaT_read_pair = ','.join(metaT_reads_list)
        #make_bowtie2_idx(mapping_ref, mapping_result_dir, threads)
        run_bowtie2(mapping_ref, each_metaT_read_pair, mapping_index, mapping_result_dir, sam_name, threads)
    elif len(metaT_reads_list) / 2 >= 2 and len(metaT_reads_list) % 2 == 0:
        #make_bowtie2_idx(mapping_ref, mapping_result_dir, threads)
        for i in range(0, len(metaT_reads_list), 2):
            j = i + 1
            sam_name = Path(metaT_reads_list[i]).stem.rsplit('_', 1)[0]
            sam_names.append(sam_name)
            each_metaT_read_pair = f'{metaT_reads_list[i]},{metaT_reads_list[j]}'
            run_bowtie2(mapping_ref, each_metaT_read_pair, mapping_index, mapping_result_dir, sam_name, threads)
    else:
        sys.exit('You input reads are not in pairs, please check') 
        
    # Step 2 Filter sam file
    for sam_name in sam_names:
        input_sam_file = f'{mapping_result_dir}/{sam_name}.sam'
        convert_sam_to_sorted_bam(input_sam_file, threads)
        out_bam_file = input_sam_file.replace('.sam', '.bam', 1)
        out_sorted_bam_file = input_sam_file.replace('.sam', '.sorted.bam', 1)
        filtered_bam_file = input_sam_file.replace('.sam', '.filtered.bam', 1)
        aligned_length = 50
        filter_sorted_bam(out_sorted_bam_file, filtered_bam_file, reads_mapping_identity_cutoff, aligned_length, threads)
        os.system(f'rm {input_sam_file} {out_bam_file} {out_sorted_bam_file}')    
        
    # Step 3 Get coverage
    bam_files = ''
    bam_files_list = []
    for sam_name in sam_names:
        bam_files_list.append(f'{mapping_result_dir}/{sam_name}.filtered.bam')
    bam_files = ' '.join(bam_files_list)
        
    os.system(f'coverm contig --methods rpkm --bam-files {bam_files} --threads {threads} > {mapping_result_dir}/all_coverm_raw_result.txt')    

# Step 1 Mapping metaT (diel_cycle_metaT) to all_phage_species_rep
mapping_ref_all_phage_species_rep = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_24hrs_metaT_ref.based_on_all_phage_species_rep.fa'
mapping_index_all_phage_species_rep  = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_24hrs_metaT_ref.based_on_all_phage_species_rep.bowtie2_idx'
metaT_reads_diel_cycle_metaT = '/storage1/Reads/TYMEFLIES_reads/diel_cycle_metaT_reads/diel_cycle_metaT_reads.rRNA_removed.filtered_1.fastq.gz,/storage1/Reads/TYMEFLIES_reads/diel_cycle_metaT_reads/diel_cycle_metaT_reads.rRNA_removed.filtered_2.fastq.gz'
mapping_result_dir_diel_cycle_metaT = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/diel_cycle_metaT___all_phage_species_rep_mapping_result_dir'
mapping_metaT_reads(mapping_ref_all_phage_species_rep, mapping_index_all_phage_species_rep, metaT_reads_diel_cycle_metaT, mapping_result_dir_diel_cycle_metaT, '0.97', '10')