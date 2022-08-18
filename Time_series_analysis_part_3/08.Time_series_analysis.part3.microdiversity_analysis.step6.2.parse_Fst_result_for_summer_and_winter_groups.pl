#!/usr/bin/perl

use strict;
use warnings;

# Aim: Parse Fst result for winter group (at the 1st place) and summer group (at the 2nd place)

# Step 1 Store Fst result of each gene
my %Fst = (); # $gene => [0] $fst [1] $pi_1 [2] $pi_2 [3] $cov_1 [4] $cov_2
open IN, "ls Summer_vs_Winter_Fst_analysis/Fst_for_each_viral_gn/*.tsv |";
while (<IN>){
	chomp;
	my $file = $_;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^gene/){
			my @tmp = split (/\t/);
			my $gene = $tmp[0];
			my $fst = $tmp[2];
			my $pi_1 = $tmp[3];
			my $pi_2 = $tmp[4];
			my $cov_1 = $tmp[5];
			my $cov_2 = $tmp[6];
			if (! $fst){
				$fst = "NA";
			}
			if (! $pi_1){
				$pi_1 = "NA";
			}
			if (! $pi_2){
				$pi_2 = "NA";
			}
			if (! $cov_1){
				$cov_1 = "NA";
			}
			if (! $cov_2){
				$cov_2 = "NA";
			}			
			$Fst{$gene}[0] = $fst;
			$Fst{$gene}[1] = $pi_1;
			$Fst{$gene}[2] = $pi_2;
			$Fst{$gene}[3] = $cov_1;
			$Fst{$gene}[4] = $cov_2;
		}
	}
	close INN;
}
close IN;

# Step 2 Store pNpS for each gene
my %Gene2pNpS = (); # $gene => [0] $pNpS_1 [1] $pNpS_2
foreach my $gene (sort keys %Fst){
	$Gene2pNpS{$gene}[0] = "NA";
	$Gene2pNpS{$gene}[1] = "NA";
}

open IN, "Summer_vs_Winter_Fst_analysis/All_winter_IMG.IS/output/All_winter_IMG.IS_gene_info.tsv";
while (<IN>){
	chomp;
	if (!/^scaffold/){
		my @tmp = split (/\t/);
		my $gene = $tmp[1];
		my $pNpS_1 = $tmp[12];
		if (! $pNpS_1){
			$pNpS_1 = "NA";
		}
		
		$Gene2pNpS{$gene}[0] = $pNpS_1;
	}
}
close IN;

open IN, "Summer_vs_Winter_Fst_analysis/All_summer_IMG.IS/output/All_summer_IMG.IS_gene_info.tsv";
while (<IN>){
	chomp;
	if (!/^scaffold/){
		my @tmp = split (/\t/);
		my $gene = $tmp[1];
		my $pNpS_2 = $tmp[12];
		if (! $pNpS_2){
			$pNpS_2 = "NA";
		}
		
		$Gene2pNpS{$gene}[1] = $pNpS_2;
	}
}
close IN;

# Step 3 Filter Fst result
# The following requirements are used:
# (1) Fst >= 0.5
# (2) pi in winter group > pi in summer group
# (3) gene N/S SNV ratio in winter group < gene N/S SNV ratio in summer group
# (4) gene coverages in winter group and summer group both > 5×

my %Fst_filtered = (); # $gene => [0] $fst [1] $pi_1 [2] $pi_2 [3] $cov_1 [4] $cov_2 [5] $pNpS_1 [6] $pNpS_2
foreach my $gene (sort keys %Fst){
	my $fst = $Fst{$gene}[0];
	my $pi_1 = $Fst{$gene}[1];
	my $pi_2 = $Fst{$gene}[2];
	my $cov_1 = $Fst{$gene}[3];
	my $cov_2 = $Fst{$gene}[4];
	my $pNpS_1 = $Gene2pNpS{$gene}[0];
	my $pNpS_2 = $Gene2pNpS{$gene}[1];
	
	my $logic = 0; # Record how many requirements are meet
	if ($fst ne "NA" and $fst >= 0.75){
		$logic++;
	}
	
	if ($pi_1 ne "NA" and $pi_2 ne "NA" and $pi_1 > $pi_2){
		$logic++;
	}
	
	if ($pNpS_1 ne "NA" and $pNpS_2 ne "NA" and $pNpS_1 < $pNpS_2){
		$logic++;
	}	
	
	if ($cov_1 ne "NA" and $cov_2 ne "NA" and $cov_1 > 5 and $cov_2 > 5){
		$logic++;
	}		
	
	if ($logic == 4){
		$Fst_filtered{$gene}[0] = $fst;
		$Fst_filtered{$gene}[1] = $pi_1;
		$Fst_filtered{$gene}[2] = $pi_2;
		$Fst_filtered{$gene}[3] = $cov_1;
		$Fst_filtered{$gene}[4] = $cov_2;		
		$Fst_filtered{$gene}[5] = $pNpS_1;	
		$Fst_filtered{$gene}[6] = $pNpS_2;			
	}
}

# Step 4 Write down filtered Fst result
open OUT, ">Summer_vs_Winter_Fst_analysis/Fst_for_each_viral_gn/All_filtered_Fst_result.txt";
print OUT "Gene\tFst\tPi_1\tPi_2\tCov_1\tCov_2\tpNpS_1\tpNpS_2\n";
foreach my $gene (sort keys %Fst_filtered){
	print OUT "$gene\t$Fst_filtered{$gene}[0]\t$Fst_filtered{$gene}[1]\t$Fst_filtered{$gene}[2]\t$Fst_filtered{$gene}[3]\t$Fst_filtered{$gene}[4]\t$Fst_filtered{$gene}[5]\t$Fst_filtered{$gene}[6]\n";
}
close OUT;