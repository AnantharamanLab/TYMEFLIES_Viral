#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(shuffle);

# Aim: Make the input file for refraction curve for species

# Step 1 Store species hash and make %IMGID2Species2presence
## Store viral genomes that are present in Species_level_vOTUs_cluster.txt
my %Species_map = (); # $gn_rep => $gns
my %IMGID2Species2presence = (); # $img_id => $gn_rep (represent for each species) => $presence (1)
open IN, "/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gn_rep = $tmp[0];
	my $gns = $tmp[1];
	$Species_map{$gn_rep} = $gns;
	# Store all the IMG ID for this species
	my @IMG_ID = ($gns=~/(33\d+?)\_\_/g);
	my %IMG_ID = map { $_ => 1 } @IMG_ID;
	foreach my $img_id (sort keys %IMG_ID){
		$IMGID2Species2presence{$img_id}{$gn_rep} = 1;
	}
}
close IN;

# Step 2 Start from a random sample and add sample gradually
## Step 2.1 Store all sample (metagenome) ID and 
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '^33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

my @IMGID = sort keys %IMGID;

## Step 2.2 Pick a random sample and add sample gradually to get the Sample_num2species_num hash, do it for 1 reps
for(my $j=1; $j<=9; $j++){
	my @IMGID_random = shuffle @IMGID; # Shuffle the array
	my %Sample_num2species_num = (); # $sample_num => $species_num
	my @Sample_included = (); # Store all the samples that should be included in the analysis
	for(my $i=0; $i<=$#IMGID_random; $i++){
		my $img_id = $IMGID_random[$i];
		push @Sample_included, $img_id;
		
		my $sample_num = scalar (@Sample_included);
		my $species_num = 0;
		
		foreach my $gn_rep (sort keys %Species_map){
			my $logic = 0;
			OUTTER: foreach my $img_id (@Sample_included){
				if (exists $IMGID2Species2presence{$img_id}{$gn_rep}){
					$logic = 1; last OUTTER;
				}
			}
			$species_num += $logic;
		}
		
		$Sample_num2species_num{$sample_num} = $species_num;
	}

	# Write down results
	my @Sample_num = (1..471);
	open OUT, ">Cluster_phage_genomes/Sample_num2species_num_for_refraction_curve.$j.txt";
	foreach my $i (@Sample_num){
		print OUT "$i\t$Sample_num2species_num{$i}\n";
	}
	close OUT;
}
