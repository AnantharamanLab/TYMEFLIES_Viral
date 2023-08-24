#!/usr/bin/perl

use strict;
use warnings;

# AIM: Grep specific AMG genes (both nucleotide and amino acid sequences) by KO ID from all metagenomes

# The variables
# The gene name and KO ID for grepping
my %Gene2ko_id_for_grepping = ();
$Gene2ko_id_for_grepping{"pmoC-amoC"} = "K10946";

# Step 1. Store AMG summary table
my %AMG_summary = (); # $pro => [0] $date_and_season [1] $amg_ko [2] $amg_ko_name
                      # [3] $pfam [4] $pfam_name [5] $metabolisms [6] $pathways
open IN, "AMG_analysis/AMG_summary.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $pro = $tmp[0];
	my $date_and_season = $tmp[1];
	my $amg_ko = $tmp[2];
	my $amg_ko_name = $tmp[3];
	my $pfam = $tmp[4];
	my $pfam_name = $tmp[5];
	my $metabolisms = $tmp[6];
	my $pathways = $tmp[7];
	$AMG_summary{$pro}[0] = $date_and_season;
	$AMG_summary{$pro}[1] = $amg_ko;
	$AMG_summary{$pro}[2] = $amg_ko_name;
	$AMG_summary{$pro}[3] = $pfam;
	$AMG_summary{$pro}[4] = $pfam_name;
	$AMG_summary{$pro}[5] = $metabolisms;
	$AMG_summary{$pro}[6] = $pathways;
}
close IN;

# Step 2. Grep fasta and faa from each metagenome
foreach my $gene_symbol (sort keys %Gene2ko_id_for_grepping){
	my $ko_id_targeted = $Gene2ko_id_for_grepping{$gene_symbol};
	
	my %Ffn = (); # Store the targeted AMG fasta for this gene
	my %Faa = (); # Store the targeted AMG faa for this gene
	
	my %IMG_ID = (); # Store the $img_id that are involved with these AMGs
	my %Targeted_pro_id = (); # Store the $pro that are targeted AMGs
	my %Bin = (); # Store the $bin (bin name) that are involved with these AMGs
	
	foreach my $pro (sort keys %AMG_summary){
		if ($AMG_summary{$pro}[1] eq $ko_id_targeted){
			my ($img_id) = $pro =~ /^(.+?)\_\_/;
			$IMG_ID{$img_id} = 1;
			my ($bin) = $pro =~ /^(.+?\_\_.+?)\_\_/;
			$Bin{$bin} = 1; # Examples: 3300046832__vRhyme_unbinned1249 or 3300046832__vRhyme_488
			$Targeted_pro_id{$pro} = 1;
		}
	}
	
	foreach my $img_id (sort keys %IMG_ID){
		# Store all the fasta files into %Ffn
		open IN, "find /storage1/data11/TYMEFLIES_phage/$img_id/vRhyme_best_bins_fasta_parsed -name '*.ffn' |";
		while (<IN>){
			chomp;
			my $file = $_;
			# Filter bin names, only store sequences in targeted bin 
			my ($bin_name) = $file =~ /vRhyme_best_bins_fasta_parsed\/(.+?)\.ffn$/; 
			if (exists $Bin{$bin_name}){
				my %Seq_tmp = _store_seq("$file");
				%Ffn = (%Seq_tmp, %Ffn);
			}
		}
		close IN;	
		
		# Store all the faa files into %Faa
		open IN, "find /storage1/data11/TYMEFLIES_phage/$img_id/vRhyme_best_bins_fasta_parsed -name '*.faa' |";
		while (<IN>){
			chomp;
			my $file = $_;
			# Filter bin names, only store sequences in targeted bin 
			my ($bin_name) = $file =~ /vRhyme_best_bins_fasta_parsed\/(.+?)\.faa$/;
			if (exists $Bin{$bin_name}){
				my %Seq_tmp = _store_seq("$file");
				%Faa = (%Seq_tmp, %Faa);
			}
		}
		close IN;		
	}
	
	# Delete non-targeted AMG hits
	foreach my $key (sort keys %Ffn){
		my ($key_new) = $key =~ /^>(.+?)$/; # The clean pro id without ">"
		if (!exists $Targeted_pro_id{$key_new}){
			delete $Ffn{$key};
			delete $Faa{$key}; # %Ffn and %Faa have the same set of keys
		}
	}
	
	# Store the Ffn and Faa files
	`mkdir Targeted_AMGs`;
	
	open OUT, ">Targeted_AMGs/$gene_symbol.$ko_id_targeted.fasta";
	foreach my $key (sort keys %Ffn){
		print OUT "$key\n$Ffn{$key}\n";
	}
	close OUT;
	
	open OUT, ">Targeted_AMGs/$gene_symbol.$ko_id_targeted.faa";
	foreach my $key (sort keys %Faa){
		print OUT "$key\n$Faa{$key}\n";
	}
	close OUT;	
}



# Subroutines
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