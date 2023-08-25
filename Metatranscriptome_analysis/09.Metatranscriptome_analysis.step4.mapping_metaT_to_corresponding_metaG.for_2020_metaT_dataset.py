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
    
# Aim: Mapping metaT to corresponding metaG for 2020 metaT datasets
# The metaT ID to corresponding metaG ID:
# ME_2020_07_24_5m_C => 3300044631 (2018/7/24)
# ME_2020_08_05_10M_B => 3300044608 (2018/8/6)
# ME_2020_08_25_10M_B => 3300042355 (2018/8/25)
# ME_2020_10_19_5M_B => 3300034116 (2018/10/24)

# Use 3300044631 (2018/7/24) as the mapping reference for ME_2020_07_24_5m_C 
# Use 3300044608 (2018/8/6) as the mapping reference for ME_2020_08_05_10M_B 
# Use 3300042355 (2018/8/25) as the mapping reference for ME_2020_08_25_10M_B 
# Use 3300034116 (2018/10/24) as the mapping reference for ME_2020_10_19_5M_B

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
 
def run_bowtie2(fasta, input_read_pair, working_dir, sam_name, num_threads):
    file_name = Path(fasta).stem
    index_name = file_name + ".bowtie2_idx"
    
    # Mapping 
    mapping_cmd = f'bowtie2 -x {working_dir}/{index_name} -1 {input_read_pair.split(",")[0]} -2 {input_read_pair.split(",")[1]} -S {working_dir}/{sam_name}.sam -p {int(num_threads)} --no-unal --quiet --mm 1> /dev/null'
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
    
def mapping_metaT_reads(mapping_ref, metaT_reads, mapping_result_dir, reads_mapping_identity_cutoff, threads):
    threads = int(threads)
    
    # Step 1 Run Bowtie2
    os.mkdir(mapping_result_dir)
    metaT_reads_list = metaT_reads.split(',')
    sam_names = []
    if len(metaT_reads_list) / 2 == 1:
        sam_name = Path(metaT_reads_list[0]).stem.rsplit('_', 1)[0]
        sam_names.append(sam_name)
        each_metaT_read_pair = ','.join(metaT_reads_list)
        make_bowtie2_idx(mapping_ref, mapping_result_dir, threads)
        run_bowtie2(mapping_ref, each_metaT_read_pair, mapping_result_dir, sam_name, threads)
    elif len(metaT_reads_list) / 2 >= 2 and len(metaT_reads_list) % 2 == 0:
        make_bowtie2_idx(mapping_ref, mapping_result_dir, threads)
        for i in range(0, len(metaT_reads_list), 2):
            j = i + 1
            sam_name = Path(metaT_reads_list[i]).stem.rsplit('_', 1)[0]
            sam_names.append(sam_name)
            each_metaT_read_pair = f'{metaT_reads_list[i]},{metaT_reads_list[j]}'
            run_bowtie2(mapping_ref, each_metaT_read_pair, mapping_result_dir, sam_name, threads)
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

# Step 1 Mapping metaT (ME_2020_07_24_5m_C) to corresponding metaG for 3300044631 
mapping_ref_3300044631 = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_2020_metaT_ref_for_3300044631.based_on_MetaG.fa'
metaT_reads_ME_2020_07_24_5m_C = '/slowdata/Reads/Paired_Bact_Vir_Mendota/Metatranscriptomes/Metatranscriptomes_processed/ME_2020_07_24_5m_C_R1.fastq,/slowdata/Reads/Paired_Bact_Vir_Mendota/Metatranscriptomes/Metatranscriptomes_processed/ME_2020_07_24_5m_C_R2.fastq'
mapping_result_dir_ME_2020_07_24_5m_C = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_07_24_5m_C___3300044631_mapping_result_dir'
mapping_metaT_reads(mapping_ref_3300044631, metaT_reads_ME_2020_07_24_5m_C, mapping_result_dir_ME_2020_07_24_5m_C, '0.97', '10')


# Step 2 Mapping metaT (ME_2020_08_05_10M_B) to corresponding metaG for 3300044608 
mapping_ref_3300044608 = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_2020_metaT_ref_for_3300044608.based_on_MetaG.fa'
metaT_reads_ME_2020_08_05_10M_B = '/slowdata/Reads/Paired_Bact_Vir_Mendota/Metatranscriptomes/Metatranscriptomes_processed/ME_2020_08_05_10M_B_R1.fastq,/slowdata/Reads/Paired_Bact_Vir_Mendota/Metatranscriptomes/Metatranscriptomes_processed/ME_2020_08_05_10M_B_R2.fastq'
mapping_result_dir_ME_2020_08_05_10M_B = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_08_05_10M_B___3300044608_mapping_result_dir'
mapping_metaT_reads(mapping_ref_3300044608, metaT_reads_ME_2020_08_05_10M_B, mapping_result_dir_ME_2020_08_05_10M_B, '0.97', '10')


# Step 3 Mapping metaT (ME_2020_08_25_10M_B) to corresponding metaG for 3300042355 
mapping_ref_3300042355 = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_2020_metaT_ref_for_3300042355.based_on_MetaG.fa'
metaT_reads_ME_2020_08_25_10M_B = '/slowdata/Reads/Paired_Bact_Vir_Mendota/Metatranscriptomes/Metatranscriptomes_processed/ME_2020_08_25_10M_B_R1.fastq,/slowdata/Reads/Paired_Bact_Vir_Mendota/Metatranscriptomes/Metatranscriptomes_processed/ME_2020_08_25_10M_B_R2.fastq'
mapping_result_dir_ME_2020_08_25_10M_B = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_08_25_10M_B___3300042355_mapping_result_dir'
mapping_metaT_reads(mapping_ref_3300042355, metaT_reads_ME_2020_08_25_10M_B, mapping_result_dir_ME_2020_08_25_10M_B, '0.97', '10')


# Step 4 Mapping metaT (ME_2020_10_19_5M_B) to corresponding metaG for 3300034116 
mapping_ref_3300034116 = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_2020_metaT_ref_for_3300034116.based_on_MetaG.fa'
metaT_reads_ME_2020_10_19_5M_B = '/slowdata/Reads/Paired_Bact_Vir_Mendota/Metatranscriptomes/Metatranscriptomes_processed/ME_2020_10_19_5M_B_R1.fastq,/slowdata/Reads/Paired_Bact_Vir_Mendota/Metatranscriptomes/Metatranscriptomes_processed/ME_2020_10_19_5M_B_R2.fastq'
mapping_result_dir_ME_2020_10_19_5M_B = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_10_19_5M_B___3300034116_mapping_result_dir'
mapping_metaT_reads(mapping_ref_3300034116, metaT_reads_ME_2020_10_19_5M_B, mapping_result_dir_ME_2020_10_19_5M_B, '0.97', '10')





