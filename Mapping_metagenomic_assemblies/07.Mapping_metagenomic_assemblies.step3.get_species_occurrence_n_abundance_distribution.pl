#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get the species-level vOTU occurrence and abundance distribution pattern

# Step 1 Store all the viral genome abundance (normalized by 100M reads/metagenome)
## Step 1.1 Store scaffold coverage for each phage genome
my %Scf2cov = (); # $scf => $cov (within individual metagenome); 
                  # here $scf did not contain the viral genome in the front
				  # like "Ga0453728_0000141_fragment_2"
                
open IN, "ls 33**/vRhyme_result/vRhyme_coverage_files/vRhyme_coverage_values.tsv | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($img_id) = $file =~ /(33\d+?)\//;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^scaffold/){
			my @tmp = split (/\t/);
			my $scf = $tmp[0];
			my $cov = $tmp[1];
			$Scf2cov{$scf} = $cov;
		}
	}
	close INN;
}
close IN;

## Step 1.2 Store viral genome to viral scaffold hash
my %Viral_gn2viral_scf = (); # $viral_gn => $viral_scfs (collection of $viral_scf, separated by "\t")
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/All_phage_genomes_headers.txt";
while (<IN>){
	chomp;
	my ($viral_scf) = $_ =~ /^>(.+?)$/;
	my ($viral_gn) = $viral_scf =~ /^(.+?\_\_.+?)\_\_/;
	if (!exists $Viral_gn2viral_scf{$viral_gn}){
		$Viral_gn2viral_scf{$viral_gn} = $viral_scf;
	}else{
		$Viral_gn2viral_scf{$viral_gn} .= "\t".$viral_scf;
	}
}
close IN;

## Step 1.3 Get read numbers for each metagenome
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

## Step 1.4 Get viral genome abundance (normalized)
my %Viral_gn2IMG2abun = (); # $viral_gn => $img_id => $abun (normalized)
foreach my $viral_gn (sort keys %Viral_gn2viral_scf){
	my ($img_id) = $viral_gn =~ /^(.+?)\_\_/;
	my $abun = 0; # the normalized viral genome abundance 
	
	my @Viral_scf_abun = (); # Store all the depth values of viral scfs
	my @Viral_scfs =  split (/\t/, $Viral_gn2viral_scf{$viral_gn});
	foreach my $viral_scf (@Viral_scfs){
		my ($viral_scf_short) = $viral_scf =~ /^.+?\_\_.+?\_\_(.+?)$/;
		my $depth = $Scf2cov{$viral_scf_short};
		push @Viral_scf_abun, $depth;
	}
	
	$abun = _avg(@Viral_scf_abun);
	
	my $read_num = $IMG_ID2read_num{$img_id};
	
	$abun = $abun / ($read_num / 100000000);
	
	$Viral_gn2IMG2abun{$viral_gn}{$img_id} = $abun;
}

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
