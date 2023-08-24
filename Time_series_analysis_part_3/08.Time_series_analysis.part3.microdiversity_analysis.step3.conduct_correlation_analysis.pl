#!/usr/bin/perl

use strict;
use warnings;

# Aim: Conduct Spearman's correlation test between species cov with nt diversity and SNP density

# Note: This script should be run under env vRhyme (conda activate vRhyme)

# Step 1 Conduct Spearman's correlation test for species containing four AMGs
## Step 1.1 Conduct Spearman's correlation test between cov and nt diversity
my %Viral_gn_four_AMG2cov_array = (); # $viral_gn => $cov_array (include the gn name as the first element)
open IN, "MetaPop/Viral_gn2IMG2cov_norm_filtered.four_AMGs.txt";
while (<IN>){
	chomp;
	if (!/^Head/){
		my @tmp = split (/\t/);
		shift @tmp; # Delete the first element
		my $viral_gn = $tmp[0];
		my $cov_array = join("\t", @tmp);
		$Viral_gn_four_AMG2cov_array{$viral_gn} = $cov_array;
	}
}
close IN;

my %Viral_gn_four_AMG2nt_diversity_array = (); # $viral_gn => $nt_diversity_array (include the gn name as the first element)
open IN, "MetaPop/Viral_gn2IMG2nt_diversity.four_AMGs.txt";
while (<IN>){
	chomp;
	if (!/^Head/){
		my @tmp = split (/\t/);
		shift @tmp; # Delete the first element
		my $viral_gn = $tmp[0];
		my $nt_diversity_array = join("\t", @tmp);
		$Viral_gn_four_AMG2nt_diversity_array{$viral_gn} = $nt_diversity_array;
	}
}
close IN;

### Write down tmp result to do the correlation test
`mkdir MetaPop/tmp_folder_for_correlation_test`;

my %Corr_result_four_AMG = (); # $viral_gn => [0] $corr_coeff [1] $pvalue
foreach my $viral_gn (sort keys %Viral_gn_four_AMG2cov_array){
	open OUT, ">MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.input.txt";
	print OUT "$Viral_gn_four_AMG2cov_array{$viral_gn}\n";
	print OUT "$Viral_gn_four_AMG2nt_diversity_array{$viral_gn}\n";
	close OUT;
	
	`python3 calc_spearman_correlation.py -i MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.input.txt -o MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.output.txt`;

	open IN, "MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.output.txt";
	while (<IN>){
		chomp;
		my @tmp = split (/\t/);
		my $pvalue = $tmp[2];
		my $corr_coeff = $tmp[3];
		$Corr_result_four_AMG{$viral_gn}[0] = $corr_coeff;
		$Corr_result_four_AMG{$viral_gn}[1] = $pvalue;
	}
	close IN;
}

`rm -r MetaPop/tmp_folder_for_correlation_test`;

## Step 1.2 Conduct Spearman's correlation test between cov and SNP density
my %Viral_gn_four_AMG2snp_density_array = (); # $viral_gn => $snp_density_array (include the gn name as the first element)
open IN, "MetaPop/Viral_gn2IMG2snp_density.four_AMGs.txt";
while (<IN>){
	chomp;
	if (!/^Head/){
		my @tmp = split (/\t/);
		shift @tmp; # Delete the first element
		my $viral_gn = $tmp[0];
		my $snp_density_array = join("\t", @tmp);
		$Viral_gn_four_AMG2snp_density_array{$viral_gn} = $snp_density_array;
	}
}
close IN;

### Write down tmp result to do the correlation test
`mkdir MetaPop/tmp_folder_for_correlation_test`;

### Store corr_coeff and pvalue for SNP diversity
# $viral_gn => [2] $corr_coeff [3] $pvalue
foreach my $viral_gn (sort keys %Viral_gn_four_AMG2cov_array){
	open OUT, ">MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.input.txt";
	print OUT "$Viral_gn_four_AMG2cov_array{$viral_gn}\n";
	print OUT "$Viral_gn_four_AMG2snp_density_array{$viral_gn}\n";
	close OUT;
	
	`python3 calc_spearman_correlation.py -i MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.input.txt -o MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.output.txt`;

	open IN, "MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.output.txt";
	while (<IN>){
		chomp;
		my @tmp = split (/\t/);
		my $pvalue = $tmp[2];
		my $corr_coeff = $tmp[3];
		$Corr_result_four_AMG{$viral_gn}[2] = $corr_coeff;
		$Corr_result_four_AMG{$viral_gn}[3] = $pvalue;
	}
	close IN;
}

`rm -r MetaPop/tmp_folder_for_correlation_test`;

## Step 1.3 Write down the correlation test result 
### Store the viral species containing four AMGs
my %Viral_species_containing_four_AMGs = (); # $viral_gn => $info (the corresponding information for each viral species gn)
my @Viral_species_containing_four_AMGs = (); # Store the $viral_gn order
open IN, "viral_species_containing_four_AMGs.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $info = $tmp[0];
	my $viral_gn = $tmp[1];
	$Viral_species_containing_four_AMGs{$viral_gn} = $info;
	
	push @Viral_species_containing_four_AMGs, $viral_gn;
}
close IN;

open OUT, ">MetaPop/Correlation_between_cov_and_nt_diversity_and_SNP_density.for_four_AMGs.txt";
print OUT "Viral gn character\tViral gn\tcorr coeff nt diversity\tp value nt diversity\tcorr coeff SNP diversity\tp value SNP diversity\n";
foreach my $viral_gn (@Viral_species_containing_four_AMGs){
	print OUT "$Viral_species_containing_four_AMGs{$viral_gn}\t$viral_gn\t$Corr_result_four_AMG{$viral_gn}[0]\t$Corr_result_four_AMG{$viral_gn}[1]\t$Corr_result_four_AMG{$viral_gn}[2]\t$Corr_result_four_AMG{$viral_gn}[3]\n";
}
close OUT;

# Step 2 Conduct Spearman's correlation test for all species 
## Step 2.1 Conduct Spearman's correlation test between cov and nt diversity
my %Viral_gn2cov_array = (); # $viral_gn => $cov_array (include the gn name as the first element)
open IN, "MetaPop/Viral_gn2IMG2cov_norm_filtered.txt";
while (<IN>){
	chomp;
	if (!/^Head/){
		my @tmp = split (/\t/);
		my $viral_gn = $tmp[0];
		my $cov_array = join("\t", @tmp);
		$Viral_gn2cov_array{$viral_gn} = $cov_array;
	}
}
close IN;

my %Viral_gn2nt_diversity_array = (); # $viral_gn => $nt_diversity_array (include the gn name as the first element)
open IN, "MetaPop/Viral_gn2IMG2nt_diversity.txt";
while (<IN>){
	chomp;
	if (!/^Head/){
		my @tmp = split (/\t/);
		my $viral_gn = $tmp[0];
		my $nt_diversity_array = join("\t", @tmp);
		$Viral_gn2nt_diversity_array{$viral_gn} = $nt_diversity_array;
	}
}
close IN;

### Write down tmp result to do the correlation test
`mkdir MetaPop/tmp_folder_for_correlation_test`;

my %Corr_result = (); # $viral_gn => [0] $corr_coeff [1] $pvalue
foreach my $viral_gn (sort keys %Viral_gn2cov_array){
	open OUT, ">MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.input.txt";
	print OUT "$Viral_gn2cov_array{$viral_gn}\n";
	print OUT "$Viral_gn2nt_diversity_array{$viral_gn}\n";
	close OUT;
	
	`python3 calc_spearman_correlation.py -i MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.input.txt -o MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.output.txt`;

	open IN, "MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.output.txt";
	while (<IN>){
		chomp;
		my @tmp = split (/\t/);
		my $pvalue = $tmp[2];
		my $corr_coeff = $tmp[3];
		$Corr_result{$viral_gn}[0] = $corr_coeff;
		$Corr_result{$viral_gn}[1] = $pvalue;
	}
	close IN;
}

`rm -r MetaPop/tmp_folder_for_correlation_test`;

## Step 2.2 Conduct Spearman's correlation test between cov and SNP density
my %Viral_gn2snp_density_array = (); # $viral_gn => $snp_density_array (include the gn name as the first element)
open IN, "MetaPop/Viral_gn2IMG2snp_density.txt";
while (<IN>){
	chomp;
	if (!/^Head/){
		my @tmp = split (/\t/);
		my $viral_gn = $tmp[0];
		my $snp_density_array = join("\t", @tmp);
		$Viral_gn2snp_density_array{$viral_gn} = $snp_density_array;
	}
}
close IN;

### Write down tmp result to do the correlation test
`mkdir MetaPop/tmp_folder_for_correlation_test`;

### Store corr_coeff and pvalue for SNP diversity
foreach my $viral_gn (sort keys %Viral_gn2cov_array){
	open OUT, ">MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.input.txt";
	print OUT "$Viral_gn2cov_array{$viral_gn}\n";
	print OUT "$Viral_gn2snp_density_array{$viral_gn}\n";
	close OUT;
	
	`python3 calc_spearman_correlation.py -i MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.input.txt -o MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.output.txt`;

	open IN, "MetaPop/tmp_folder_for_correlation_test/$viral_gn.corr_test.output.txt";
	while (<IN>){
		chomp;
		my @tmp = split (/\t/);
		my $pvalue = $tmp[2];
		my $corr_coeff = $tmp[3];
		$Corr_result{$viral_gn}[2] = $corr_coeff;
		$Corr_result{$viral_gn}[3] = $pvalue;
	}
	close IN;
}

`rm -r MetaPop/tmp_folder_for_correlation_test`;

## Step 2.3 Write down the correlation test result 
open OUT, ">MetaPop/Correlation_between_cov_and_nt_diversity_and_SNP_density.txt";
print OUT "Viral gn\tcorr coeff nt diversity\tp value nt diversity\tcorr coeff SNP diversity\tp value SNP diversity\n";
foreach my $viral_gn (sort keys %Viral_gn2cov_array){
	print OUT "$viral_gn\t$Corr_result{$viral_gn}[0]\t$Corr_result{$viral_gn}[1]\t$Corr_result{$viral_gn}[2]\t$Corr_result{$viral_gn}[3]\n";
}
close OUT;