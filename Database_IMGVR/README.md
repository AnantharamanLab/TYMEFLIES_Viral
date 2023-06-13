# IMG/VR v4.1

IMG/VR database v4.1 release Dec. 2022  (for [Cluster vMAGs](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Cluster_vMAGs))

The link: https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html

- NAME: IMGVR_all_nucleotides.fna
It contains 15,722,824 sequences

- NAME: IMGVR_all_Host_information.tsv (The host information of a subset of viral genomes)
It contains 15,677,623 lines. For each genome, there might be multiple host information by multiple host prediction methods

- NAME: IMGVR_all_Sequence_information.tsv (All viral genome metadata)
It contains 15,677,623 phage genomes, and they are clustered into 8,606,551 vOTUs (dereplicated)

- NAME: IMGVR_all_phage_vOTU_representatives.txt (The vOTU representative map)
  It contains 3,216,119 vOTUs and 3,216,119 vOTU representatives. Representatives are picked firstly by genome quality and then by genome size

  

[script] 01.grep_individual_phage_genomes.py

Grep sequences in the "IMGVR_all_nucleotides.fna" and divide them into individual phage genomes.

Step 1. Check phage genome name congruency.

Step 2. Write down "IMGVR_all_genome_scaffolds.txt" to store scaffolds to genome map

[output] IMGVR_all_genome_scaffolds.txt.gz



[script] 02.make_vOTU_representatives_list.py

Make vOTU representative list; pick the high quality and long phage genome from each vOTU as the representative genome; store individual OTU representative sequences

[output] IMGVR_all_phage_vOTU_representatives.txt.gz









