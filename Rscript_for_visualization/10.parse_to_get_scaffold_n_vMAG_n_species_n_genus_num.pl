#!/usr/bin/perl

use strict;
use warnings;

# AIM: Get the viral scaffold number
#              vMAG number
#			   species-level vOTU number
#			   genus-level vOTU number

# Step 1 Get the viral scaffold number
my $viral_scaffold_num = 0;

open IN, "cat /storage1/data11/TYMEFLIES_phage/33*/VIBRANT_33*.a.v2.min2000/phage_results.min2000.txt |";
while (<IN>){
	chomp;
	if (/^Ga/){
		$viral_scaffold_num++;
	}
}
close IN;

# Step 2 Get the vMAG number
my $vMAG_num = 0;

open IN, "cat /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed/Each_bin_info.txt |";
while (<IN>){
	chomp;
	if (/^33/){
		$vMAG_num++;
	}
}
close IN;

# Step 3 Get the species-level vOTU number

my $species_level_vOTU_num = 0;

open IN, "/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt";
while (<IN>){
	chomp;
	if (/^33/){
		$species_level_vOTU_num++;
	}
}
close IN;

# Step 4 Get the genus-level vOTU number

my $genus_level_vOTU_num = 0;

## Add genus from genus_clusters.txt
my %All_virus_in_genus_cluster = (); # $gn => 1
open IN, "/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/genus_clusters.txt";
while (<IN>){
	chomp;
	if (/^33/){
		$genus_level_vOTU_num++;
		my @tmp = split (/\t/);
		foreach my $key (@tmp){
			$All_virus_in_genus_cluster{$key} = 1;
		}
	}
}
close IN;

## Add genus that is not from genus_clusters.txt
open IN, "cat /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed/Each_bin_info.txt |";
while (<IN>){
	chomp;
	if (/^33/){
		my @tmp = split (/\t/);
		my $gn = $tmp[0];
		if (! exists $All_virus_in_genus_cluster{$gn}){
			$genus_level_vOTU_num++;
		}
	}
}
close IN;

# Step 5 Write down result
open OUT, ">Scaffold_vMAG_Species_Genus_num.txt";
print OUT "Viral scaffold number\t$viral_scaffold_num\n";
print OUT "vMAG number\t$vMAG_num\n";
print OUT "Species-level vOTU number\t$species_level_vOTU_num\n";
print OUT "Genus-level vOTU number\t$genus_level_vOTU_num\n";
close OUT;

