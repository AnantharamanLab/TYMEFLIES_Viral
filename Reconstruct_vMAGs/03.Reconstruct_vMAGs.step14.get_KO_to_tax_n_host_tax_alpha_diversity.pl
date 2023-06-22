#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(shuffle);

# AIM: Get the KO to tax and host tax alpha diversity input tables

# Step 1 Store all KO summary information
my %Pro2ko = (); # $pro => $ko
my %KO2pros = (); # $ko => $pros (collection of $pro separated by "\t")
my %KO2viral_gns = (); # $ko => $viral_gns (collection of $viral_gn separated by "\t")
open IN, "AMG_analysis/AMG_summary.txt";
while (<IN>){
	chomp;
	if (!/^Protein/){
		my @tmp = split (/\t/);
		my $pro = $tmp[0];
		my $ko = $tmp[2];
		$Pro2ko{$pro} = $ko;
		if (!exists $KO2pros{$ko}){
			$KO2pros{$ko} = $pro;
		}else{
			$KO2pros{$ko} .= "\,".$pro;
		}
		
		my ($viral_gn) = $pro =~ /^(.+?\_\_.+?)\_\_/;
		if (!exists $KO2viral_gns{$ko}){
			$KO2viral_gns{$ko} = $viral_gn;
		}else{
			$KO2viral_gns{$ko} .= "\t".$viral_gn;
		}
	}
}
close IN;

# Step 2 Get KO to viral gn number hash
my %KO2viral_gn_num = (); # $ko => $viral_gn_num
foreach my $ko (sort keys %KO2viral_gns){
	my $viral_gn_num = 0;
	my @Viral_gns = split (/\t/, $KO2viral_gns{$ko});
	$viral_gn_num = scalar @Viral_gns;
	$KO2viral_gn_num{$ko} = $viral_gn_num;
}
close OUT;

# Step 3 Get the alpha diversity input table for viral tax
## Step 3.1 Store the viral gn to tax hash
my %Viral_gn2tax = (); # $gn => $tax;
open IN, "/storage1/data11/TYMEFLIES_phage/Taxonomic_classification/Each_bin_tax_combined_result.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$Viral_gn2tax{$tmp[0]} = $tmp[1];
}
close IN;

## Step 3.2 Store the KO2family2viral_gn_num hash
my %KO2family2viral_gn_num = (); # $ko => $family (also contains order in the front) => $viral_gn_num
my %KO2viral_gns_w_informative_family_num = (); # $ko => $viral_gns_w_informative_family_num
my %Family_for_viral_tax = (); # $family => 1
my %KO2viral_gns_w_informative_family = (); # $ko => $viral_gns_w_informative_family
foreach my $ko (sort keys %KO2viral_gn_num){
	my @Viral_gns_w_informative_family = (); # The total viral genomes for informative family (excluding "Unclassified", "NA;NA")
	
	my @Viral_gns = split (/\t/, $KO2viral_gns{$ko});
	foreach my $viral_gn (@Viral_gns){
		if (exists $Viral_gn2tax{$viral_gn}){
			my $tax = $Viral_gn2tax{$viral_gn};
			my @tmp = split (/\;/, $tax);
			my $family = $tmp[4]."\;".$tmp[5]; # Only store the order and family
			if ($family ne "NA;NA"){
				push @Viral_gns_w_informative_family, $viral_gn;
				$KO2family2viral_gn_num{$ko}{$family}++;
				$Family_for_viral_tax{$family} = 1;
			}
		}
	}
	
	my $viral_gns_w_informative_family_num = scalar @Viral_gns_w_informative_family;
	$KO2viral_gns_w_informative_family_num{$ko} = $viral_gns_w_informative_family_num;
	$KO2viral_gns_w_informative_family{$ko} = join("\t", @Viral_gns_w_informative_family);
}

## Step 3.3 Store the KO2family2viral_gn_num2 hash, only randomly takes 100 viral genomes with informative family assigned
my %KO2family2viral_gn_num2 = (); # $ko => $family (also contains order in the front) => $viral_gn_num
my %Family_for_viral_tax2 = (); # $family => 1
foreach my $ko (sort keys %KO2viral_gns_w_informative_family_num){
	if ($KO2viral_gns_w_informative_family_num{$ko} >= 100){
		my @Viral_gns_w_informative_family = split (/\t/,$KO2viral_gns_w_informative_family{$ko});
		@Viral_gns_w_informative_family = shuffle @Viral_gns_w_informative_family;
		splice @Viral_gns_w_informative_family, 100; # Get the first 100 elements
		
		foreach my $viral_gn (@Viral_gns_w_informative_family){
			my $tax = $Viral_gn2tax{$viral_gn};
			my @tmp = split (/\;/, $tax);
			my $family = $tmp[4]."\;".$tmp[5]; # Only store the order and family
			$Family_for_viral_tax2{$family} = 1;
			$KO2family2viral_gn_num2{$ko}{$family}++;
		}
	}
}

## Step 3.4 Write down KO2family2viral_gn_num2 hash
open OUT, ">KO2family2viral_gn_num.txt";
my $row=join("\t", sort keys %Family_for_viral_tax2);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %KO2family2viral_gn_num2){
			print OUT $tmp1."\t";
			my @tmp = ();
			foreach my $tmp2 (sort keys %Family_for_viral_tax2) {       
					if (exists $KO2family2viral_gn_num2{$tmp1}{$tmp2}){
							push @tmp, $KO2family2viral_gn_num2{$tmp1}{$tmp2};
					}else{              
							push @tmp,"0";
					}
			}
			print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 4 Get the alpha diversity input table for viral host tax
## Step 4.1 Store the viral gn to host tax hash
my %Viral_gn2host_tax = (); # $gn => $tax;
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/Viral_gn2host_tax_final.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$Viral_gn2host_tax{$tmp[0]} = $tmp[1];
}
close IN;

## Step 4.2 Store the KO2host_family2viral_gn_num hash
my %KO2host_family2viral_gn_num = (); # $ko => $host_family (also contains order in the front) => $viral_gn_num
my %KO2viral_gns_w_informative_host_family_num = (); # $ko => $viral_gns_w_informative_host_family_num
my %Host_family_for_viral_tax = (); # $host_family => 1
my %KO2viral_gns_w_informative_host_family = (); # $ko => $viral_gns_w_informative_host_family
foreach my $ko (sort keys %KO2viral_gn_num){
	my @Viral_gns_w_informative_host_family = (); # The total viral genomes for informative host family (excluding "Unclassified", "o__;f__")
	
	my @Viral_gns = split (/\t/, $KO2viral_gns{$ko});
	foreach my $viral_gn (@Viral_gns){
		if (exists $Viral_gn2host_tax{$viral_gn}){
			my $host_tax = $Viral_gn2host_tax{$viral_gn};
			my @tmp = split (/\;/, $host_tax);
			my $host_family = $tmp[3]."\;".$tmp[4]; # Only store the order and family of host
			if ($host_family ne "o__;f__"){
				push @Viral_gns_w_informative_host_family, $viral_gn;
				$KO2host_family2viral_gn_num{$ko}{$host_family}++;
				$Host_family_for_viral_tax{$host_family} = 1;
			}
		}
	}
	
	my $viral_gns_w_informative_host_family_num = scalar @Viral_gns_w_informative_host_family;
	$KO2viral_gns_w_informative_host_family_num{$ko} = $viral_gns_w_informative_host_family_num;
	$KO2viral_gns_w_informative_host_family{$ko} = join("\t", @Viral_gns_w_informative_host_family);
}

## Step 4.3 Store the KO2host_family2viral_gn_num2 hash, only randomly takes 25 viral genomes with informative host family assigned
my %KO2host_family2viral_gn_num2 = (); # $ko => $host_family (also contains order in the front) => $viral_gn_num
my %Host_family_for_viral_tax2 = (); # $host_family => 1
foreach my $ko (sort keys %KO2viral_gns_w_informative_host_family_num){
	if ($KO2viral_gns_w_informative_host_family_num{$ko} >= 25){
		my @Viral_gns_w_informative_host_family = split (/\t/,$KO2viral_gns_w_informative_host_family{$ko});
		@Viral_gns_w_informative_host_family = shuffle @Viral_gns_w_informative_host_family;
		splice @Viral_gns_w_informative_host_family, 25; # Get the first 25 elements
		
		foreach my $viral_gn (@Viral_gns_w_informative_host_family){
			my $host_tax = $Viral_gn2host_tax{$viral_gn};
			my @tmp = split (/\;/, $host_tax);
			my $host_family = $tmp[3]."\;".$tmp[4]; # Only store the order and family of host
			$Host_family_for_viral_tax2{$host_family} = 1;
			$KO2host_family2viral_gn_num2{$ko}{$host_family}++;
		}
	}
}

## Step 4.4 Write down KO2host_family2viral_gn_num2 hash
open OUT, ">KO2host_family2viral_gn_num.txt";
my $row2=join("\t", sort keys %Host_family_for_viral_tax2);
print OUT "Head\t$row2\n";
foreach my $tmp1 (sort keys %KO2host_family2viral_gn_num2){
			print OUT $tmp1."\t";
			my @tmp = ();
			foreach my $tmp2 (sort keys %Host_family_for_viral_tax2) {       
					if (exists $KO2host_family2viral_gn_num2{$tmp1}{$tmp2}){
							push @tmp, $KO2host_family2viral_gn_num2{$tmp1}{$tmp2};
					}else{              
							push @tmp,"0";
					}
			}
			print OUT join("\t",@tmp)."\n";
}
close OUT;
