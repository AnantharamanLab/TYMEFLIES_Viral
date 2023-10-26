#!/usr/bin/perl

use strict;
use warnings;

# Aim: Run metapop for all "*.viral_species_rep_containing_AMG.id90.bam" files (the viral species representative genomes containing AMG)
# This script should be run under the conda env - metapop (use "conda activate metapop" to activate the env)

`metapop -i ./viral_species_rep_bams.for_summer_vs_winter -r ./reference_fasta_for_metapop -g All_phage_species_rep_gn_containing_AMG.mdfed.genes --id_min 93 --threads 20 --snp_scale both --no_macro --no_viz  -o MetaPop.for_summer_vs_winter > metapop.for_summer_vs_winter.log`;

# id_min set to 93% (the species genome boundary as suggested)
# Use genes provided by me

