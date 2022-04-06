#!/usr/bin/perl

use strict;
use warnings;

# Aim: Parse metapop result to get AMG coverage
# This script should be run under the conda env - "conda activate /home/kieft/miniconda3/envs/pyscripts"

# Step 1 Get and write down the AMG coordinates
## Step 1.1 Store AMG summary table 
my %AMG_summary = (); # $pro => $ko
my %KOs= (); # $ko => 1;
my %IMG2date = (); # $img_id => $date_n_season
open IN, "AMG_analysis/AMG_summary.txt";
while (<IN>){
	chomp;
	if (!/^Pro/){
		my @tmp = split (/\t/);
		my $pro = $tmp[0];
		my $ko = $tmp[2];
		my $ko_detail = $tmp[3];
		my $date_n_season = $tmp[1];
		$AMG_summary{$pro} = $ko;
		my ($img_id) = $pro =~ /^(33.+?)\_/;
		$IMG2date{$img_id} = $date_n_season;
		$KOs{$ko} = $ko_detail;
	}
}
close IN;

## Step 1.2 Store and write down AMG coordinates 
my %AMG_coordinates = (); # $amg => $scf, $amg, $start, $stop connected by "\t"
open IN, "All_phage_species_rep_gn_containing_AMG.genes";
while (<IN>){
	chomp;
	if (/^>/){
		my $line = $_;
		my ($amg, $start, $stop) = $line =~ /^>(.+?) \# (.+?) \# (.+?) \#/;
		my ($scf) = $amg =~ /^(.+)\_/;
		if (exists $AMG_summary{$amg}){
			$AMG_coordinates{$amg} = "$scf\t$amg\t$start\t$stop";
		}
	}
}
close IN;

open OUT, ">MetaPop/All_phage_species_rep_gn_containing_AMG_coordinates.txt";
print OUT "scaffold\tregion\tstart\tstop\n";
foreach my $amg (sort keys %AMG_coordinates){
	print OUT "$AMG_coordinates{$amg}\n";
}
close OUT;

## Step 1.3 Make prophage coordinates file
`cat /storage1/data11/TYMEFLIES_phage/33*/VIBRANT_33*.a/VIBRANT_results_33*.a/VIBRANT_integrated_prophage_coordinates_33*.a.tsv > MetaPop/VIBRANT_prophage_coordinates.txt`;

# Step 2 Write down batch run command and run it
`mkdir MetaPop/AMG_coverage_result`;
open OUT, ">tmp.batch_run_to_get_AMG_coverage.sh";
open IN, "ls /storage1/data11/TYMEFLIES_phage/MetaPop/04.Depth_per_Pos/*.tsv |";
while (<IN>){
	chomp;
	my $depth_file = $_;
	my ($basename) = $depth_file =~ /(33*.+?\.viral_species_rep\.id90)/;
	print OUT "/slowdata/scripts/python_scripts/cov_by_region.py -i $depth_file -r MetaPop/All_phage_species_rep_gn_containing_AMG_coordinates.txt -c MetaPop/VIBRANT_prophage_coordinates.txt -f MetaPop/01.Genomes_and_Genes/all_genomes.fasta -o MetaPop/AMG_coverage_result/$basename.AMG_cov.txt --no_header\n";
}
close IN;
close OUT;

`cat tmp.batch_run_to_get_AMG_coverage.sh | parallel -j 20`;

`rm tmp.batch_run_to_get_AMG_coverage.sh`;
