#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get the species-level vOTU occurrence and abundance distribution pattern


# Step 1 Store all the viral genome abundance (normalized by 100M reads/metagenome)
## Step 1.1 Get read numbers for each metagenome
my %IMG_ID2read_num = ();
open IN, "TYMEFLIES_metagenome_info.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $img_id = $tmp[0];
		my $read_num = $tmp[13];
		$IMG_ID2read_num{$img_id} = $read_num;
	}
}	
close IN;

## Step 1.2 Store the previously calculated viral_gn abun (normalized)
## The viral_gn abun (normalized) file "viral_gn2depth_normalized.txt" was 
## calculated by "03.Reconstruct_vMAGs.step7.get_all_virus_abundance.py"              
my %Viral_gn2IMG2abun = (); # $viral_gn => $img_id => $abun (normalized)
open IN, "viral_gn2depth_normalized.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $viral_gn = $tmp[0];
	my $abun_normalized = $tmp[1];
	my ($img_id) = $viral_gn =~ /^(\d+?)_/;
	$Viral_gn2IMG2abun{$viral_gn}{$img_id} = $abun_normalized;
}	
close IN;

# Step 2 Get the species abundance hash
## Step 2.1 Store the species hash
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

## Step 2.2 Store %Species2IMG2abun hash
my %Species2IMG2abun = (); # $species => $img_id => $abun
foreach my $species (sort keys %Species){
	foreach my $img_id (sort keys %IMG_ID2read_num){
		my $abun = 0;
		
		my @Gns = split (/\,/, $Species{$species});
		foreach my $gn (@Gns){
			my ($img_id2) = $gn =~ /^(.+?)\_\_/;
			if ($img_id eq $img_id2){
				if ($Viral_gn2IMG2abun{$gn}{$img_id}){
					$abun += $Viral_gn2IMG2abun{$gn}{$img_id};
				}
			}
		}
		
		$Species2IMG2abun{$species}{$img_id} = $abun;
	}
}

# Step 3 Store species occurrence
my %Species2IMG = (); # $species => $img_ids (collection of $img_id, separated by "\t")
foreach my $gn_rep (sort keys %Species){
	my @Gns = split (/\,/, $Species{$gn_rep});
	
	foreach my $gn (@Gns){
		my ($img_id) = $gn =~ /^(.+?)\_\_/;
		if (!exists $Species2IMG{$gn_rep}){
			$Species2IMG{$gn_rep} = $img_id;
		}else{
			$Species2IMG{$gn_rep} .= "\t".$img_id;
		}
	}
}

# Step 4 Store and write down %Species2occurrence_n_abundance hash 
my %Species2occurrence_n_abundance = (); # $species => [0] $occurrence [1] $abundance
foreach my $species (sort keys %Species2IMG){
	my $occurrence = 0;
	my @IMG = split (/\t/, $Species2IMG{$species});
	$occurrence = scalar @IMG;
	
	my $abundance = 0; # Store the mean abundance in IMG that species occur
	my @Abundance = (); # Store abundance in IMG that species occur
	foreach my $img_id (@IMG){
		my $abun = $Species2IMG2abun{$species}{$img_id};
		push @Abundance, $abun;
	}
	$abundance = _avg(@Abundance);
	
	$Species2occurrence_n_abundance{$species}[0] = $occurrence;
	$Species2occurrence_n_abundance{$species}[1] = $abundance;
}

open OUT, ">Species2occurrence_n_abundance.txt";
foreach my $species (sort keys %Species2occurrence_n_abundance){
	print OUT "$species\t$Species2occurrence_n_abundance{$species}[0]\t$Species2occurrence_n_abundance{$species}[1]\n";
}
close OUT;



# Subroutine 
sub _avg {
    my $total;
    $total += $_ foreach @_;
    # sum divided by number of components.
    return $total / @_;
}
