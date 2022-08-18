#!/usr/bin/perl

use strict;
use warnings;

# Aim: Parse metapop result to get AMG coverage for each year
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

### Change the old gene to new gene
my %Old_gene2new_gene_map = (); # $gene_old => $gene_new
open IN, "New_gene2old_gene_map.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gene_new = $tmp[0]; my $gene_old = $tmp[1];
	$Old_gene2new_gene_map{$gene_old} = $gene_new;
}
close IN;

foreach my $pro (sort keys %AMG_summary){
	if ($Old_gene2new_gene_map{$pro}){
		my $ko = $AMG_summary{$pro};
		my $gene_new = $Old_gene2new_gene_map{$pro};
		delete $AMG_summary{$pro}; # Delete the old gene and its value
		$AMG_summary{$gene_new} = $ko; # Add the new gene and its value
	}
}

## Step 1.2 Store and write down AMG coordinates 
my %AMG_coordinates = (); # $amg => $scf, $amg, $start, $stop connected by "\t"
open IN, "All_phage_species_rep_gn_containing_AMG.mdfed.genes"; # Note: Using mdfed gene files here 
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

open OUT, ">MetaPop.for_each_year/MetaPop/All_phage_species_rep_gn_containing_AMG_coordinates.txt";
print OUT "scaffold\tregion\tstart\tstop\n";
foreach my $amg (sort keys %AMG_coordinates){
	print OUT "$AMG_coordinates{$amg}\n";
}
close OUT;

=pod
## Do not need prophage coordinates now, due to that prophage coordinates have been changed already - all start positions are changed to 1, and start protein IDs are changed to 1
## Step 1.3 Make prophage coordinates file
`cat /storage1/data11/TYMEFLIES_phage/33*/VIBRANT_33*.a/VIBRANT_results_33*.a/VIBRANT_integrated_prophage_coordinates_33*.a.tsv > MetaPop/VIBRANT_prophage_coordinates.txt`;
=cut

# Step 2 Write down batch run command and run it
`mkdir MetaPop.for_each_year/MetaPop/AMG_coverage_result`;
open OUT, ">tmp.batch_run_to_get_AMG_coverage.sh";
open IN, "ls /storage1/data11/TYMEFLIES_phage/MetaPop.for_each_year/MetaPop/04.Depth_per_Pos/*.tsv |";
while (<IN>){
	chomp;
	my $depth_file = $_;
	my ($basename) = $depth_file =~ /(20*.+?\.viral_species_rep\.id90)/;
	print OUT "/storage1/data14/for_chao/cov_by_region.py -i $depth_file -r MetaPop.for_each_year/MetaPop/All_phage_species_rep_gn_containing_AMG_coordinates.txt -f MetaPop.for_each_year/MetaPop/01.Genomes_and_Genes/all_genomes.fasta -o MetaPop.for_each_year/MetaPop/AMG_coverage_result/$basename.AMG_cov.txt --no_header\n";
}
close IN;
close OUT;

`mkdir tmp_dir`;

`cat tmp.batch_run_to_get_AMG_coverage.sh | parallel -j 30 --tmpdir tmp_dir`;

`rm tmp.batch_run_to_get_AMG_coverage.sh`;

`rm -r tmp_dir`;