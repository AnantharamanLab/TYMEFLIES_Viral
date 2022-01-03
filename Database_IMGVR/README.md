# Database processing

**1 IMG V/R v3**

IMG_VR_2020-10-12_5.1 - IMG/VR v3  (for [Cluster vMAGs](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Cluster_vMAGs))

The link: https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html

- NAME: IMGVR_all_nucleotides.fna  (The orignal nt sequences of IMG_VR_2020-10-12_5.1 - IMG/VR v3)
It contains 2377994 sequences

- NAME: IMGVR_all_nucleotides_header.txt (The head line (without ">") of IMGVR_all_nucleotides.fna)

- NAME: IMGVR_all_Host_information.tsv (The host information of a subset of viral genomes)
It contains 2053759 lines. For each genome, there might be multiple host information by multiple host prediction methods

- NAME: IMGVR_all_Sequence_information.tsv (All viral genome metadata)
It contains 2332702 phage genomes, and they are clustered into 935362 vOTUs (dereplicated)

- NAME: IMGVR_all_genome_scaffolds.txt

  It contains three columns: "Phage genome ID", "Number of scaffolds in the genome", and "Scaffold ID" (separated by "\t")

- NAME: IMGVR_all_phage_vOTU_representatives.txt (The vOTU representative map)
It contains 935362 vOTUs and 935362 vOTU representatives. Representatives are picked firstly by genome quality and then by genome size



[script] 01.grep_individual_phage_genomes.pl

Grep sequences in the "IMGVR_all_nucleotides.fna" and divide them into individual phage genomes.

Step 1. Check phage genome name congruency.

Step 2. Make individual phage genomes in a new folder.

[output] IMGVR_all_genome_scaffolds.txt.gz



[script] 02.make_vOTU_representatives_list.pl

Make vOTU representative list; pick the high quality and long phage genome from each vOTU as the representative genome; store individual OTU representative sequences

[output] IMGVR_all_phage_vOTU_representatives.txt.gz



2 NCBI RefSeq viral

1) NCBI RefSeq viral: Downloaded from https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/

The last modified date of this database is: 2021-11-04 (for Taxonomical Classification)

2) ICTV_Master_Species_List_2020.v1.txt is the 2021/5/10 release version (MSL #36)





The "viral.protein.ictv_8_rank_tax.txt" is the tax map for RefSeq viral proteins (not all proteins have a corresponding tax tag)
in which the tax is in the style of ICTV (8 ranks in total)









