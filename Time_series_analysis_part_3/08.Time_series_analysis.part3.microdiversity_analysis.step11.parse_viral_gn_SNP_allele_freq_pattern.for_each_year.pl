#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum);

# Aim: Parse the viral species (species rep gn) SNP allele frequencies for each year
#      Get the linear regression patterns for each viral gn
#      and the Spearman' rank correlation result for each viral gn

# Step 1 Store %All_SNP_pos2year2allele_freq hash
my %All_SNP_pos2year2allele_freq = (); # $snp_pos => $year => $allele_freq
my @Head = ();
my %Year = (); # $year => 1
my %All_SNP_pos2viral_gn = (); # $snp_pos => $viral_gn
open IN, "MetaPop.for_each_year/MetaPop/All_SNP_pos2year2allele_freq.txt";
while (<IN>){
	chomp;
	if (/^Head/){
		my @tmp = split(/\t/);
		@Head = @tmp;
	}else{
		my @tmp = split(/\t/);
		for(my $i=1; $i<=$#tmp; $i++){
			my $snp_pos = $tmp[0];
			my $allele_freq = $tmp[$i];
			my $year = $Head[$i];
			my ($viral_gn) = $snp_pos =~ /^(.+?)\_\_Ga/;
			$Year{$year} = 1;
			$All_SNP_pos2year2allele_freq{$snp_pos}{$year} = $allele_freq;
			$All_SNP_pos2viral_gn{$snp_pos} = $viral_gn;
		}
	}
}
close IN;

my %Viral_gn2snp_pos = (); # $viral_gn => collection of $snp_pos
foreach my $snp_pos (sort keys %All_SNP_pos2viral_gn){
	my $viral_gn = $All_SNP_pos2viral_gn{$snp_pos};
	if (!exists $Viral_gn2snp_pos{$viral_gn}){
		$Viral_gn2snp_pos{$viral_gn} = $snp_pos;
	}else{
		$Viral_gn2snp_pos{$viral_gn} .= "\t".$snp_pos;
	}
}

# Step 2 Store %Viral_gn2year2allele_freq_mean hash
my %Viral_gn2year2allele_freq_mean = (); # $viral_gn => $year => $allele_freq_mean
foreach my $viral_gn (sort keys %Viral_gn2snp_pos){
	foreach my $year (sort keys %Year){
		my @SNP_pos = split (/\t/, $Viral_gn2snp_pos{$viral_gn});
		
		my $allele_freq_mean = "NA";
		my @Allele_freq = (); # Store all $allele_freq within @SNP_pos
		foreach my $snp_pos (@SNP_pos){
			my $allele_freq = $All_SNP_pos2year2allele_freq{$snp_pos}{$year};
			if ($allele_freq ne "NA"){
				push @Allele_freq, $allele_freq;
			}
		}
		if (@Allele_freq){
			$allele_freq_mean = _mean(@Allele_freq);
		}
		
		$Viral_gn2year2allele_freq_mean{$viral_gn}{$year} = $allele_freq_mean;
	}
}

# Step 3 Write down %Viral_gn2year2allele_freq_mean hash
open OUT, ">MetaPop.for_each_year/MetaPop/Viral_gn2year2allele_freq_mean.txt";
my $row=join("\t", sort keys %Year);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %Viral_gn2snp_pos){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Year) {       
                if (exists $Viral_gn2year2allele_freq_mean{$tmp1}{$tmp2}){
                        push @tmp, $Viral_gn2year2allele_freq_mean{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"NA";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 4 Get linear regression result and Spearman' rank correlation result
`mkdir MetaPop.for_each_year/MetaPop/Viral_gn_allele_freq_tmp_dir`;

## Step 4.1 Get linear regression result
open OUTT, ">MetaPop.for_each_year/MetaPop/batch_input_for_allele_freq_pattern.txt";
my %Viral_gn2regression_parameters = (); # $viral_gn => $regression_parameters
foreach my $viral_gn (sort keys %Viral_gn2year2allele_freq_mean){
	my $line_viral_gn_and_allele_freq = "";
	my @Line_viral_gn_and_allele_freq = ();
	push @Line_viral_gn_and_allele_freq, $viral_gn;
	foreach my $year (sort keys %Year){
		push @Line_viral_gn_and_allele_freq, $Viral_gn2year2allele_freq_mean{$viral_gn}{$year};
	}
	$line_viral_gn_and_allele_freq = join("\t", @Line_viral_gn_and_allele_freq);
	open OUT, ">MetaPop.for_each_year/MetaPop/Viral_gn_allele_freq_tmp_dir/$viral_gn.allele_freq_input.txt";
	print OUT $line_viral_gn_and_allele_freq."\n";
	close OUT;
	print OUTT "MetaPop.for_each_year/MetaPop/Viral_gn_allele_freq_tmp_dir/$viral_gn.allele_freq_input.txt\tMetaPop.for_each_year/MetaPop/Viral_gn_allele_freq_tmp_dir/$viral_gn.allele_freq_output.txt\n";
}
close OUTT;

`python3 /storage1/data14/for_chao/calc_regression.py -i MetaPop.for_each_year/MetaPop/batch_input_for_allele_freq_pattern.txt --batch-mode`;

open IN, "find MetaPop.for_each_year/MetaPop/Viral_gn_allele_freq_tmp_dir/ -name '*.allele_freq_output.txt' | ";
while (<IN>){
	chomp;
	my $file = $_;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (/^33/){
			my @tmp = split (/\t/);
			my $viral_gn = $tmp[0];
			my $reg_slope = $tmp[1];
			my $reg_yint = $tmp[2];
			my $reg_rsquared = "NA";
			if ($tmp[3]){
				$reg_rsquared = $tmp[3];
			}
		
			$Viral_gn2regression_parameters{$viral_gn} = "$reg_slope\t$reg_yint\t$reg_rsquared";
		}
	}
	close INN;
}
close IN;

## Step 4.2 Get Spearman' rank correlation result
open OUTT, ">MetaPop.for_each_year/MetaPop/batch_input_for_allele_freq_pattern.for_spearman_test.txt";
my %Viral_gn2spearman_parameters = (); # $viral_gn => $spearman_parameters
foreach my $viral_gn (sort keys %Viral_gn2year2allele_freq_mean){
	my $line_viral_gn_and_allele_freq = "";
	my @Line_viral_gn_and_allele_freq = ();
	push @Line_viral_gn_and_allele_freq, $viral_gn;
	foreach my $year (sort keys %Year){
		push @Line_viral_gn_and_allele_freq, $Viral_gn2year2allele_freq_mean{$viral_gn}{$year};
	}
	$line_viral_gn_and_allele_freq = join("\t", @Line_viral_gn_and_allele_freq);
	open OUT, ">MetaPop.for_each_year/MetaPop/Viral_gn_allele_freq_tmp_dir/$viral_gn.allele_freq_input.for_spearman_test.txt";
	print OUT "Head\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\t11\t12\t13\t14\t15\t16\t17\t18\t19\t20\n";
	print OUT $line_viral_gn_and_allele_freq."\n";
	close OUT;
	print OUTT "MetaPop.for_each_year/MetaPop/Viral_gn_allele_freq_tmp_dir/$viral_gn.allele_freq_input.for_spearman_test.txt\tMetaPop.for_each_year/MetaPop/Viral_gn_allele_freq_tmp_dir/$viral_gn.allele_freq_output.for_spearman_test.txt\n";
}
close OUTT;

`python3 /storage1/data14/for_chao/calc_spearman_correlation.py -i MetaPop.for_each_year/MetaPop/batch_input_for_allele_freq_pattern.for_spearman_test.txt --batch-mode`;

open IN, "find MetaPop.for_each_year/MetaPop/Viral_gn_allele_freq_tmp_dir/ -name '*.allele_freq_output.for_spearman_test.txt' | ";
while (<IN>){
	chomp;
	my $file = $_;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (/^Head/){
			my @tmp = split (/\t/);
			my $viral_gn = $tmp[1];
			my $spearman_pval = $tmp[2];
			my $spearman_corr = $tmp[3];
			$Viral_gn2spearman_parameters{$viral_gn} = "$spearman_corr\t$spearman_pval";
		}
	}
	close INN;
}
close IN;

`rm -r MetaPop.for_each_year/MetaPop/Viral_gn_allele_freq_tmp_dir`;

# Step 4.3 Write down linear regression and Spearman' rank correlation result
open OUT, ">MetaPop.for_each_year/MetaPop/Viral_gn2regression_and_spearman_parameters.txt";
foreach my $viral_gn (sort keys %Viral_gn2regression_parameters){
	print OUT "$viral_gn\t$Viral_gn2regression_parameters{$viral_gn}\t$Viral_gn2spearman_parameters{$viral_gn}\n";
}
close OUT;

## Store the viral species containing four AMGs
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

open OUT, ">MetaPop.for_each_year/MetaPop/Viral_gn2regression_and_spearman_parameters.for_four_AMGs.txt";
foreach my $viral_gn (sort keys %Viral_gn2regression_parameters){
	if (exists $Viral_species_containing_four_AMGs{$viral_gn}){
		print OUT "$viral_gn\t$Viral_gn2regression_parameters{$viral_gn}\t$Viral_gn2spearman_parameters{$viral_gn}\n";
	}
}
close OUT;



# Subroutine
sub _mean {
    return sum(@_)/@_;
}
