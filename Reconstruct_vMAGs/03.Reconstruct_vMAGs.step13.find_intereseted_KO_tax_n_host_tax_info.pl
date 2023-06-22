#!/usr/bin/env perl

use strict;
use warnings;
use Array::Split qw( split_into );

# AIM: Find interested KO tax and host tax info (the abundance fraction)

# Step 1 Store all KO summary information
my %Pro2ko = (); # $pro => $ko
my %KO2pro = (); # $ko => $pros (collection of $pro separated by "\,")
open IN, "AMG_analysis/AMG_summary.txt";
while (<IN>){
	chomp;
	if (!/^Protein/){
		my @tmp = split (/\t/);
		my $pro = $tmp[0];
		my $ko = $tmp[2];
		$Pro2ko{$pro} = $ko;
		if (!exists $KO2pro{$ko}){
			$KO2pro{$ko} = $pro;
		}else{
			$KO2pro{$ko} .= "\,".$pro;
		}
	}
}
close IN;

# Step 2 Store the abundance (normalized by 100M reads/metagenome) of each protein
## Step 2.1 Get read numbers for each metagenome
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

## Step 2.2 Store scaffold coverage 
my %Scf2cov = (); # $scf => $cov (within individual metagenome)
open IN, "ls /storage1/data11/TYMEFLIES_phage/*/*.id97.coverm_depth.txt | ";
while (<IN>){
	chomp;
	my $file = $_;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^contigName/){
			my @tmp = split (/\t/);
			my $scf = $tmp[0];
			my $cov = $tmp[2];
			$Scf2cov{$scf} = $cov;
		}
	}
	close INN;
}
close IN;

## Step 2.3 Store prophage scaffold coverage (containing 'fragment' in scaffold header)
my %Prophage_scf2cov = (); # $scf => $cov (within individual metagenome)
open IN, "ls /storage1/data11/TYMEFLIES_phage/*/PropagAtE_result/PropagAtE_result.tsv | ";
while (<IN>){
	chomp;
	my $file = $_;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^prophage/){
			my @tmp = split (/\t/);
			my $prophage_scf = $tmp[0];
			my $cov = $tmp[7];
			if ($cov eq 'NA'){
				$cov = 0;
			}
			$Prophage_scf2cov{$prophage_scf} = $cov;
		}
	}
	close INN;
}
close IN;

## Step 2.4 Store pro abundance (normalized by 100M reads/metagenome)
my %Pro2abun = (); # $pro => $abun (normalized)
foreach my $pro (sort keys %Pro2ko){
	my ($img_id,$scf) = $pro =~ /^(33\d+?)\_\_vRhyme\_.+?\_\_(Ga.+)\_\d+?$/;
	my $abun = 0;
	my $abun_normalized = 0;
	if ($scf =~ /fragment/){ # If this scaffold is from prophage 
		print "$scf\n"; # Print out the $scf to make a check
		$abun = $Prophage_scf2cov{$scf}; 
	}else{		
		$abun = $Scf2cov{$scf}; # If this scaffold is not from prophage 
	}
	my $read_num = $IMG_ID2read_num{$img_id};
	if ($abun){
		$abun_normalized = $abun / ($read_num / 100000000);
	}
	$Pro2abun{$pro} = $abun_normalized;	
}

# Step 3 Get the abundance fraction of tax for each KO
## Step 3.1 Store the viral gn to tax hash
my %Viral_gn2tax = (); # $gn => $tax;
open IN, "/storage1/data11/TYMEFLIES_phage/Taxonomic_classification/Each_bin_tax_combined_result.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$Viral_gn2tax{$tmp[0]} = $tmp[1];
}
close IN;

## Step 3.2 Make the hash
my %KO2tax2abun_fraction = (); # $ko => $tax (family level) => $abun_fraction
my %Tax = (); # Store all the tax
foreach my $ko (sort keys %KO2pro){
	my @Pros = split (/\,/, $KO2pro{$ko});
	foreach my $pro (@Pros){
		my ($viral_gn) = $pro =~ /^(.+?\_\_.+?)\_\_/;
		my $tax = ""; 
		if (!exists $Viral_gn2tax{$viral_gn}){
			$tax = "Unclassified";
			$Tax{$tax} = 1;
			$KO2tax2abun_fraction{$ko}{$tax} += $Pro2abun{$pro};
		}else{
			$tax = $Viral_gn2tax{$viral_gn};
			my @tmp = split (/\;/, $tax);
			$tax = $tmp[4]."\;".$tmp[5]; # Only store the order and family
			$Tax{$tax} = 1;			
			$KO2tax2abun_fraction{$ko}{$tax} += $Pro2abun{$pro};
		}
	}
}

## Step 3.3 Write down the result
open OUT, ">AMG_analysis/KO2tax2abun_fraction.txt";
my $row=join("\t", sort keys %KO2tax2abun_fraction);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %Tax){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %KO2tax2abun_fraction) {       
                if (exists $KO2tax2abun_fraction{$tmp2}{$tmp1}){
                        push @tmp, $KO2tax2abun_fraction{$tmp2}{$tmp1};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 4 Get the abundance fraction of host tax for each KO
## Step 4.1 Store the viral gn to host tax hash
my %Viral_gn2host_tax = (); # $gn => $tax;
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/Viral_gn2host_tax_final.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$Viral_gn2host_tax{$tmp[0]} = $tmp[1];
}
close IN;

## Step 4.2 Make the hash
my %KO2host_tax2abun_fraction = (); # $ko => $host_tax (family level) => $abun_fraction
my %Host_tax = (); # Store all the host tax
foreach my $ko (sort keys %KO2pro){
	my @Pros = split (/\,/, $KO2pro{$ko});
	foreach my $pro (@Pros){
		my ($viral_gn) = $pro =~ /^(.+?\_\_.+?)\_\_/;
		my $host_tax = ""; 
		if (!exists $Viral_gn2host_tax{$viral_gn}){
			$host_tax = "Unclassified";
			$Host_tax{$host_tax} = 1;
			$KO2host_tax2abun_fraction{$ko}{$host_tax} += $Pro2abun{$pro};
		}else{
			$host_tax = $Viral_gn2host_tax{$viral_gn};
			my @tmp = split (/\;/, $host_tax);
			$host_tax = $tmp[3]."\;".$tmp[4]; # Only store the order and family 
			$Host_tax{$host_tax} = 1;			
			$KO2host_tax2abun_fraction{$ko}{$host_tax} += $Pro2abun{$pro};
		}
	}
}

## Step 4.3 Write down the result
open OUT, ">AMG_analysis/KO2host_tax2abun_fraction.txt";
my $row2=join("\t", sort keys %KO2host_tax2abun_fraction);
print OUT "Head\t$row2\n";
foreach my $tmp1 (sort keys %Host_tax){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %KO2host_tax2abun_fraction) {       
                if (exists $KO2host_tax2abun_fraction{$tmp2}{$tmp1}){
                        push @tmp, $KO2host_tax2abun_fraction{$tmp2}{$tmp1};
                }else{              
                        push @tmp,"0";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;
