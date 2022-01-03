# Taxonomic classification

(The method was based on *Nucleic Acids Res*. 2021 Jan 8;49(D1):D764-D775. doi: 10.1093/nar/gkaa946. Link: https://academic.oup.com/nar/article/49/D1/D764/5952208)

**1 Run diamond to NCBI RefSeq to determine phage genome classification**

Use all phage faa files to compare against NCBI Viral RefSeq proteins with diamond.

To see whether >= 30% of the proteins for a faa (phage genome) have a hit to Viral RefSeq, only count phage genomes that are above this cutoff.

Get the consensus taxonomic affiliation based on the best hits of individual proteins (>= 50 majority rule) for a phage genome.

[script] 15.taxonomic_classification.run_diamond_to_NCBI_RefSeq_viral.pl

**2 Run hmmsearch to VOG marker HMM to determine phage genome classification**

