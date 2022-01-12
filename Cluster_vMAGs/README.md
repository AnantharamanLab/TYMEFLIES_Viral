# Cluster vMAGs

**1 Cluster vMAGs**

Using dRep to cluster phage genomes in this study 

This script conducts phage genome clustering by two steps:

1) Cluster phage genomes in each year (20 years in total)

2) Pick vOTU representatives from each year and re-cluster by dRep to get final vOTUs

Phage genomes were placed together as the input genomes for "dRep dereplicate".

For "dRep dereplicate" settings for clustering phage genomes, please see the annotation information within the script.

[script] 04.Cluster_vMAGs.step1.run_dRep_to_make_vOTUs.pl

