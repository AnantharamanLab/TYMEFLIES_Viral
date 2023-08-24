#!/usr/bin/env python

from Bio import SeqIO
import pysam
from multiprocessing import Pool

# Aim: Filter the reads in "/storage1/Reads/TYMEFLIES_reads/diel_cycle_metaT_reads/*.fastq" that are mapped on pFN18A_DNA_transcript  

def filter_fastq(args):
    bam_file_path1, bam_file_path2, scaffold, fastq1_path, fastq2_path, output_prefix = args

    bam_file1 = pysam.AlignmentFile(bam_file_path1, "rb")
    bam_file2 = pysam.AlignmentFile(bam_file_path2, "rb")

    # Get a list of read names that were mapped to the specified scaffold in both bam files
    read_names1 = [read.query_name for read in bam_file1.fetch(scaffold)]
    read_names2 = [read.query_name for read in bam_file2.fetch(scaffold)]
    read_names = set(read_names1 + read_names2)

    # Open the original fastq files
    fastq1 = SeqIO.parse(fastq1_path, "fastq")
    fastq2 = SeqIO.parse(fastq2_path, "fastq")

    # Write the filtered reads to new fastq files
    filtered_fastq1 = open(f"{output_prefix}_1.fastq", "w")
    filtered_fastq2 = open(f"{output_prefix}_2.fastq", "w")

    # Loop through the fastq files and write the reads that do not match the list of read names
    for record1, record2 in zip(fastq1, fastq2):
        read_name = record1.id.split(" ")[0].split("/")[0][1:]
        if read_name not in read_names:
            SeqIO.write(record1, filtered_fastq1, "fastq")
            SeqIO.write(record2, filtered_fastq2, "fastq")

    # Close the fastq files
    filtered_fastq1.close()
    filtered_fastq2.close()

def main():
    bam_file_path1 = "/storage1/data11/TYMEFLIES_phage/Metatranscriptome/diel_cycle_metaT___24hrs_metaT_ref_mapping_result_dir/diel_cycle_metaT_reads_1.rRNA.filtered.bam"
    bam_file_path2 = "/storage1/data11/TYMEFLIES_phage/Metatranscriptome/diel_cycle_metaT___all_phage_species_rep_mapping_result_dir/diel_cycle_metaT_reads_1.rRNA.filtered.bam"
    scaffold = "pFN18A_DNA_transcript"
    
    # Create index files for the bam files
    pysam.index(bam_file_path1)
    pysam.index(bam_file_path2)
    
    fastq1_path = "/storage1/Reads/TYMEFLIES_reads/diel_cycle_metaT_reads/diel_cycle_metaT_reads_1.rRNA_removed.fastq"
    fastq2_path = "/storage1/Reads/TYMEFLIES_reads/diel_cycle_metaT_reads/diel_cycle_metaT_reads_2.rRNA_removed.fastq"
    output_prefix = "/storage1/Reads/TYMEFLIES_reads/diel_cycle_metaT_reads/diel_cycle_metaT_reads.rRNA_removed.filtered"

    with Pool(10) as p:
        p.map(filter_fastq, [(bam_file_path1, bam_file_path2, scaffold, fastq1_path, fastq2_path, output_prefix)])

if __name__ == '__main__':
    main()


