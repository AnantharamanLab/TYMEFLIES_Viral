# GEM

GEM  (A genomic catalog of Earth’s microbiomes) 2020-11 release (for [Host prediction](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Host_prediction)) 

A catalog of 52,515 metagenome-assembled genomes from over 10,450 metagenomes collected from diverse microbiomes to capture extant microbial metabolic and functional potential. MAGs from the GEMs catalog all meet the medium-quality level of the MIMAG standard (mean completeness = 83%, mean contamination = 1.3%), and include 9,143 assigned as high-quality based on the presence of a near-full complement of rRNAs, tRNAs, and single-copy protein coding genes.

The entire genome catalog and associated data is available at: https://portal.nersc.gov/GEM/

In total, there are 52515 MAGs (containing 8086855 contigs) in this database.



**1 Download the database**

[script] 01.download_files.sh

**2 Get statistics for all MAGs**

Parse "genome_metadata.tsv" to get information containing: TDB tax, Genome completeness, Genome contamination, and scaffolds (separated by ",").

[script] 02.get_statistics_for_all_MAGs.pl