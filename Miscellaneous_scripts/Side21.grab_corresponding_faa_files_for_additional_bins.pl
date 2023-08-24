#!/usr/bin/perl

use strict;
use warnings;

# Aim: Grab the corresponding faa files for additional MAGs from 49 metagenomes

# Step 1 Find the $img_id of 49 left metagenomes
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '^33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

my %IMGID_left = %IMGID; # $img_id => 1
open IN, "/storage1/data11/TYMEFLIES_phage/Binning_Data/TYMEFLIES_all_MAGs_stat.txt";
while (<IN>){
	chomp;
	if (!/IMG/){
		my @tmp = split (/\t/);
		my ($img_id) = $tmp[0] =~ /^(.+?)\_/;
		delete $IMGID_left{$img_id};
	}
}
close IN;

# Step 2 Store Bin to scaffolds information for all additional MAGs (730 MAGs from 49 metagenomes)
my %Bin2scaffolds = (); # $bin => $scaffolds (collection of $scaffold, separated by ",")
open IN, "/storage1/data11/TYMEFLIES_phage/Binning_Data/TYMEFLIES_all_MAGs_stat.additional_49_metagenomes.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $bin = $tmp[0];
	my $scaffolds = $tmp[4];
	$Bin2scaffolds{$bin} = $scaffolds;
}
close IN;

# Step 3 Grab each bin's faa sequences from its corresponding metagenome folder
foreach my $img (sort keys %IMGID_left){
	# Store metagenome faa 
	my %Metagenome_faa = _store_seq("$img/$img.a.faa");
	
	# Store metagenome scaffold to protein ID hash
	my %Metagenome_scaffold_to_protein = (); # $scaffold => $proteins (collection of $protein, separated by "\t")
	foreach my $key (sort keys %Metagenome_faa){
		my ($protein) = $key =~ /^>(.+?)$/;
		my ($scaffold) = $protein =~ /^(.+?\_.+?)\_/;
		if (! exists $Metagenome_scaffold_to_protein{$scaffold}){
			$Metagenome_scaffold_to_protein{$scaffold} = $protein;
		}else{
			$Metagenome_scaffold_to_protein{$scaffold} .= "\t".$protein;
		}
	}
	
	# For each bin within this metagenome
	foreach my $bin (sort keys %Bin2scaffolds){
		if ($bin =~ /$img\_/){ # $bin belonging to this metagenome
			my @Scaffolds = split (/\,/, $Bin2scaffolds{$bin});
			# Store all faa sequences belonging to this bin
			my %Bin_faa = ();
			foreach my $scaffold (@Scaffolds){
				my @Protein_ID = split (/\t/, $Metagenome_scaffold_to_protein{$scaffold});
				foreach my $protein (@Protein_ID){
					my $protein_w_array = ">".$protein;
					$Bin_faa{$protein_w_array} = $Metagenome_faa{$protein_w_array};
				}
			}
			
			# Write down this bin
			open OUT, ">Binning_Data/$img/$bin.faa";
			foreach my $key (sort keys %Bin_faa){
				print OUT "$key\n$Bin_faa{$key}\n";
			}
			close OUT;
		}
	}
}



# Subroutine
sub _store_seq{
	my $file = $_[0];
	my %Seq = (); my $head = "";
	open _IN, "$file";
	while (<_IN>){
		chomp;
		if (/>/){
			if (/\s/){
				($head) = $_ =~ /^(>.+?)\s/;
				$Seq{$head} = "";
			}else{
				($head) = $_ =~ /^(>.+?)$/;
				$Seq{$head} = "";
			}
		}else{
			$Seq{$head} .= $_;
		}
	}
	close _IN;
	return %Seq;
}





