#!/usr/bin/perl

use strict;
use warnings;

# Aim: Perform MCL-based clustering
# Note: run "conda activate mcl" first to set up the conda env

`mcl Cluster_phage_genomes/Shared_protein_and_AAI.genus_edges.tsv -te 8 -I 2.0 --abc -o Cluster_phage_genomes/genus_clusters.txt`;

`mcl Cluster_phage_genomes/Shared_protein_and_AAI.family_edges.tsv -te 8 -I 1.2 --abc -o Cluster_phage_genomes/family_clusters.txt`;




