#!/usr/bin/perl

use strict;
use warnings;

# Aim: Run metapop for all "*.viral_species_rep.id90.bam" files (the viral species representative genomes that contain AMGs)
# This script should be run under the conda env - metapop (use "conda activate metapop" to activate the env)

`metapop -i ./viral_species_rep_bams -r ./reference_fasta_for_metapop -g All_phage_species_rep_gn_containing_AMG.genes --norm Read_count_file_for_metapop.txt --id_min 95 --threads 30 --snp_scale both --genome_detection_cutoff 70`;

# id_min set to 95% (the species genome boundary as suggested)
# Use genes provided by me
# Percent of bases that must be covered for a sequence to be considered detected for macrodiversity analyses - 70%

