#!/usr/bin/perl

use strict;
use warnings;

# Aim: Parse to get MCL input (filter edges)
# Keep edges between genomes with >=40% AAI and genomes with either 16 shared genes or at least 20% of shared genes (relative to both genomes)
`python filter_aai.py --in_aai Cluster_phage_genomes/Shared_protein_and_AAI.txt --min_percent_shared 20 --min_num_shared 16 --min_aai 40 --out_tsv Cluster_phage_genomes/Shared_protein_and_AAI.genus_edges.tsv`;

# Keep edges between genomes with >=20% AAI and genomes with either 8 shared genes or at least 10% of shared genes (relative to both genomes)
`python filter_aai.py --in_aai Cluster_phage_genomes/Shared_protein_and_AAI.txt --min_percent_shared 10 --min_num_shared 8 --min_aai 20 --out_tsv Cluster_phage_genomes/Shared_protein_and_AAI.family_edges.tsv`;





