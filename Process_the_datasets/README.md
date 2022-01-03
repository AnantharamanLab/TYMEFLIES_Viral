# Processing the datasets

**1** **Copy and rename the fastq file according to the metagenome IMG ID**

The "Reads2IMG_ID_map.txt" contains the IMG ID, project description, and reads folder name in IMG download deposit. Use this map to rename fastq files.

[script] 01.copy_n_rename_fastq.pl

[input] Reads2IMG_ID_map.txt

**2 Get reads number and nucleotide base number for each fastq (interleaved)**

Use "reformat.sh" firstly get all the information of each interleaved fastq file, then parse it to get reads number and nucleotide base number for each fastq.

[script] 02.get_reads_num_n_base_num.pl

**3 Move IMG tar.gz and Binning Data tar.gz files**

For each IMG ID there is a IMG tar.gz file ("33***.tar.gz") contains all the assemblies and relevant annotation files inside. And for some of metagenomes, metagenome binning has been done by IMG annotation pipeline, and a "Binning_Data.tar.gz" was provided containing all the bins.

[script] 03.copy_tar.gz_files_n_dzip.pl  -  for IMG tar.gz

[script] 03.copy_tar.gz_files_n_dzip.02.pl - for Binning Data tar.gz

**4 Run bbmap to get covstat files for some of the metagenomes**

For some of metagenomes, in the folder ("33***.tar.gz") there is not a "$img_id.a\.depth\.txt" file containing the depth information of all assemblies. We run bbmap to calculate "$img_id.covstat" for these metagenomes.

[script] 04.run_bbmap_mapping.pl

