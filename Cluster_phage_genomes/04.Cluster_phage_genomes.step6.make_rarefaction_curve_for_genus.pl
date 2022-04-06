#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(shuffle);

# Aim: Make the input file for refraction curve for genus

# Step 1 Store genus hash and make %IMGID2Genus2presence
## Store viral genomes that are present in genus_cluster.txt
my %Genus_map = (); # $gn_rep => $gns
my %IMGID2Genus2presence = (); # $img_id => $gn_rep (represent for each genus) => $presence (1)
my %All_virus_in_genus_cluster = (); # $gn => 1
open IN, "/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/genus_clusters.txt";
while (<IN>){
	chomp;
	my $line = $_;
	my @tmp = split (/\t/);
	my $gn_rep = $tmp[0];
	my $gns = $line;
	$Genus_map{$gn_rep} = $gns;
	# Store all the IMG ID for this genus
	my @IMG_ID = ($gns=~/(33\d+?)\_\_/g);
	my %IMG_ID = map { $_ => 1 } @IMG_ID;
	foreach my $img_id (sort keys %IMG_ID){
		$IMGID2Genus2presence{$img_id}{$gn_rep} = 1;
	}
	
	foreach my $key (@tmp){
		$All_virus_in_genus_cluster{$key} = 1;
	}
}
close IN;

## Store viral genomes that are not present in genus_cluster.txt
open IN, "cat /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed/Each_bin_info.txt |";
while (<IN>){
	chomp;
	if (/^33/){
		my @tmp = split (/\t/);
		my $gn = $tmp[0];
		if (! exists $All_virus_in_genus_cluster{$gn}){
			$Genus_map{$gn} = $gn;
			my ($img_id) = $gn =~ /^(33.+?)\_\_/;
			$IMGID2Genus2presence{$img_id}{$gn} = 1;
		}
	}
}
close IN;


# Step 2 Start from a random sample and add sample gradually
## Step 2.1 Store all sample (metagenome) ID
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '^33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

my @IMGID = sort keys %IMGID;

## Step 2.2 Pick a random sample and add sample gradually to get the Sample_num2genus_num hash, do it for 49 reps
for(my $j=0; $j<=48; $j++){
	my @IMGID_random = shuffle @IMGID; # Shuffle the array
	my %Sample_num2genus_num = (); # $sample_num => $genus_num
	my @Sample_included = (); # Store all the samples that should be included in the analysis
	for(my $i=0; $i<=$#IMGID_random; $i++){
		my $img_id = $IMGID_random[$i];
		push @Sample_included, $img_id;
		
		my $sample_num = scalar (@Sample_included);
		my $genus_num = 0;
		
		foreach my $gn_rep (sort keys %Genus_map){
			my $logic = 0;
			OUTTER: foreach my $img_id (@Sample_included){
				if (exists $IMGID2Genus2presence{$img_id}{$gn_rep}){
					$logic = 1; last OUTTER;
				}
			}
			$genus_num += $logic;
		}
		
		$Sample_num2genus_num{$sample_num} = $genus_num;
		#print "$sample_num\t$genus_num\n";
	}

	# Write down results
	my @Sample_num = (1..465);
	open OUT, ">Cluster_phage_genomes/Sample_num2genus_num_for_refraction_curve.$j.txt";
	foreach my $i (@Sample_num){
		print OUT "$i\t$Sample_num2genus_num{$i}\n";
	}
	close OUT;
}
