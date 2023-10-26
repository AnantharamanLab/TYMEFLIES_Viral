#!/usr/bin/perl

use strict;
use warnings;

# Aim: Run metapop for all "*.viral_species_rep_containing_AMG.id90.bam" files (the viral species representative genomes containing AMG)
# This script should be run under the conda env - metapop (use "conda activate metapop" to activate the env)

# Step 1 Merge bam files for 2000-2003 and 2016-2019
`samtools merge viral_species_rep_bams.for_each_year/2000-2003.viral_species_rep_containing_AMG.id90.bam viral_species_rep_bams.for_each_year/2000.viral_species_rep_containing_AMG.id90.bam viral_species_rep_bams.for_each_year/2001.viral_species_rep_containing_AMG.id90.bam viral_species_rep_bams.for_each_year/2002.viral_species_rep_containing_AMG.id90.bam viral_species_rep_bams.for_each_year/2003.viral_species_rep_containing_AMG.id90.bam`;
`samtools merge viral_species_rep_bams.for_each_year/2016-2019.viral_species_rep_containing_AMG.id90.bam viral_species_rep_bams.for_each_year/2016.viral_species_rep_containing_AMG.id90.bam viral_species_rep_bams.for_each_year/2017.viral_species_rep_containing_AMG.id90.bam viral_species_rep_bams.for_each_year/2018.viral_species_rep_containing_AMG.id90.bam viral_species_rep_bams.for_each_year/2019.viral_species_rep_containing_AMG.id90.bam`;

## Move the merged bam files into the new folder
`mkdir viral_species_rep_bams.for_2000-2003_vs_2016-2019`;
`mv viral_species_rep_bams.for_each_year/2000-2003.viral_species_rep_containing_AMG.id90.bam viral_species_rep_bams.for_2000-2003_vs_2016-2019/2000-2003.viral_species_rep_containing_AMG.id90.bam`;
`mv viral_species_rep_bams.for_each_year/2016-2019.viral_species_rep_containing_AMG.id90.bam viral_species_rep_bams.for_2000-2003_vs_2016-2019/2016-2019.viral_species_rep_containing_AMG.id90.bam`;

# Step 2 Conduct metapop analyis
`metapop -i ./viral_species_rep_bams.for_2000-2003_vs_2016-2019 -r ./reference_fasta_for_metapop -g All_phage_species_rep_gn_containing_AMG.mdfed.genes --id_min 93 --threads 20 --snp_scale both --no_macro --no_viz  -o MetaPop.for_2000-2003_vs_2016-2019 > metapop.for_2000-2003_vs_2016-2019.log`;
# id_min set to 93% (the species genome boundary as suggested)
# Use genes provided by me

