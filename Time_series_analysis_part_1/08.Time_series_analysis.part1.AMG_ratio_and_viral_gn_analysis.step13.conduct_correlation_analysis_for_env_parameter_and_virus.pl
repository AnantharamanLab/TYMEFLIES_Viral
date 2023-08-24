#!/usr/bin/perl

use strict;
use warnings;

# Aim: Сonduct correlation analysis for environmental parameters and viruses

# Note: This script should be run under env vRhyme (conda activate vRhyme)

# Step 1 Store the environmental parameter and viruses input table
my %Env_para_and_virus = (); # $row_name => $row
my @Row_name = (); # Store all the row names
open IN, "MetaPop/Correlation_between_env_parameter_and_viruses/Env_parameter_and_viruses_input.txt";
while (<IN>){
	chomp;
	if (!/^Time/){
		my @tmp = split (/\t/);
		my $row_name = $tmp[0];
		my @tmp_wo_row_name = @tmp;
		shift @tmp_wo_row_name;
		my $row = join("\t", @tmp_wo_row_name);
		$Env_para_and_virus{$row_name} = $row;
		
		push @Row_name, $row_name;
	}
}
close IN;

# Step 2 Conduct Spearman's correlation test
my %Corr_result = (); # $row_name_1 \t $row_name_2 => [0] $corr_coeff [1] $pvalue
my $i = 1; # Store the row name pair
`mkdir MetaPop/tmp_folder_for_correlation_test`;
for(my $j=0; $j<=$#Row_name; $j++){
	for(my $k=$j+1; $k<=$#Row_name; $k++){
		my $row_name_1 = $Row_name[$j];
		my $row_name_2 = $Row_name[$k];
		
		open OUT, ">MetaPop/tmp_folder_for_correlation_test/$i.corr_test.input.txt";
		print OUT "$row_name_1\t$Env_para_and_virus{$row_name_1}\n";
		print OUT "$row_name_2\t$Env_para_and_virus{$row_name_2}\n";
		close OUT;
		
		`python3 /storage1/data14/for_chao/calc_spearman_correlation.py -i MetaPop/tmp_folder_for_correlation_test/$i.corr_test.input.txt -o MetaPop/tmp_folder_for_correlation_test/$i.corr_test.output.txt`;

		open IN, "MetaPop/tmp_folder_for_correlation_test/$i.corr_test.output.txt";
		while (<IN>){
			chomp;
			if (!/^var1/){
				my @tmp = split (/\t/);
				my $row_name_1 = $tmp[0];
				my $row_name_2 = $tmp[1];
				my $pvalue = $tmp[2];
				my $corr_coeff = $tmp[3];
				$Corr_result{"$row_name_1\t$row_name_2"}[0] = $corr_coeff;
				$Corr_result{"$row_name_1\t$row_name_2"}[1] = $pvalue;
			}
		}
		close IN;
		$i++;		
		
	}
}
#`rm -r MetaPop/tmp_folder_for_correlation_test`;

## Step 3 Write down the correlation test result 
open OUT, ">MetaPop/Correlation_between_env_parameter_and_viruses.txt";
print OUT "Item 1\tItem 2\tcorr coeff\tp value\n";
foreach my $key (sort keys %Corr_result){
	print OUT "$key\t$Corr_result{$key}[0]\t$Corr_result{$key}[1]\n";
}
close OUT;



# Subroutine 
sub _avg {
    my $total;
    $total += $_ foreach @_;
    # sum divided by number of components.
    return $total / @_;
}