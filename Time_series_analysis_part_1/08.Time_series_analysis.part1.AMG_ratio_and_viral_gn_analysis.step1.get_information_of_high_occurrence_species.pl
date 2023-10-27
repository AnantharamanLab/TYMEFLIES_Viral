#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get the information of high occurrence species (>= 40)
# Information contains: 1) Tax 2) Host tax 3) AMG KO within

# Step 1 Store the species to occurrence and abundance hash
my %Species2occurrence_n_abundance = (); # $species => [0] $occurrence [1] $abundance
open IN, "./AMG_analysis/Species2occurrence_n_abundance.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $species = $tmp[0];
	my $occurrence = $tmp[1];
	my $abundance = $tmp[2];
	$Species2occurrence_n_abundance{$species}[0] = $occurrence;
	$Species2occurrence_n_abundance{$species}[1] = $abundance;
}
close IN;

# Step 2 Store the species to all genomes within
my %Species = (); # $gn_rep => $gns (Only keep species that contain >= 10 genomes)
open IN, "/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gn_rep = $tmp[0];
	my $gns = $tmp[1];
	
	my @Gns = split (/\,/, $tmp[1]);
	if ((scalar @Gns) >= 10){
		$Species{$gn_rep} = $gns;
	}
}
close IN;

# Step 3 Store all the information for viral genomes
## Step 3.1 Store tax information
my %Viral_gn2tax = (); # $viral_gn => $tax
open IN, "/storage1/data11/TYMEFLIES_phage/Taxonomic_classification/Each_bin_tax_combined_result.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $viral_gn = $tmp[0];
	my $tax = $tmp[1];
	$Viral_gn2tax{$viral_gn} = $tax;
}
close IN;

## Step 3.2 Store host tax information
my %Viral_gn2host_tax = (); # $viral_gn => $host_tax
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/Viral_gn2host_tax_final.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $viral_gn = $tmp[0];
	my $host_tax = $tmp[1];
	$Viral_gn2host_tax{$viral_gn} = $host_tax;
}
close IN;

## Step 3.3 Store AMG KO information
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

my %Viral_gn2AMG_KOs = (); # $viral_gn => $kos (separated by "\t")
foreach my $pro (sort keys %AMG_summary){
	my ($viral_gn) = $pro =~ /^(.+?\_\_.+?)\_\_/;
	my $ko = $AMG_summary{$pro};
	if (!exists $Viral_gn2AMG_KOs{$viral_gn}){
		$Viral_gn2AMG_KOs{$viral_gn} = $ko;
	}elsif (exists $Viral_gn2AMG_KOs{$viral_gn} and $Viral_gn2AMG_KOs{$viral_gn} !~ $ko){
		$Viral_gn2AMG_KOs{$viral_gn} .= "\t".$ko;
	}
}

# Step 4 Store the information of high occurrence species
my %High_occurrence_species2info = (); # $species => [0] $tax [1] $host_tax [2] $kos 
foreach my $species (sort keys %Species2occurrence_n_abundance){
	my $occurrence = $Species2occurrence_n_abundance{$species}[0];
	if ($occurrence >= 40){
		my @Gns = split (/\,/, $Species{$species});
		
		# Get tax information
		my $tax = "";
		if (exists $Viral_gn2tax{$species}){
			$tax = $Viral_gn2tax{$species};
		}
		
		# Get host tax information
		my $host_tax = "";
		my %Host_tax = (); # Collection of $host_tax
		if (exists $Viral_gn2host_tax{$species}){
			my $host_tax_for_this_species = $Viral_gn2host_tax{$species};
			$Host_tax{$host_tax_for_this_species} = 1;
		}else{
			foreach my $gn (@Gns){
				if (exists $Viral_gn2host_tax{$gn}){
					my $host_tax_for_this_gn = $Viral_gn2host_tax{$gn};
					$Host_tax{$host_tax_for_this_gn} = 1;
				}
			}
		}
		
		if (%Host_tax){
			$host_tax = join(" \| ", sort keys (%Host_tax));
		}
		
		# Get KOs
		my $kos = "";
		my %KOs = (); # Collection of $ko
		foreach my $gn (@Gns){
			if (exists $Viral_gn2AMG_KOs{$gn}){
				my @KOs = split (/\t/, $Viral_gn2AMG_KOs{$gn});
				foreach my $ko (@KOs){
					$KOs{$ko} = 1;
				}
			}
		}
		
		if (%KOs){
			$kos = join(" \| ", sort keys (%KOs));
		}
		
		# Store all three information
		$High_occurrence_species2info{$species}[0] = $tax;
		$High_occurrence_species2info{$species}[1] = $host_tax;
		$High_occurrence_species2info{$species}[2] = $kos;
	}
}

# Step 5 Write down the information of high occurrence species
open OUT, ">High_occurrence_species2info.txt";
foreach my $species (sort keys %High_occurrence_species2info){
	my $occurrence = $Species2occurrence_n_abundance{$species}[0];
	
	my $tax = $High_occurrence_species2info{$species}[0];
	my $host_tax = $High_occurrence_species2info{$species}[1];
	my $kos = $High_occurrence_species2info{$species}[2];
	print OUT "$species\t$occurrence\t$tax\t$host_tax\t$kos\n";
}
close OUT;