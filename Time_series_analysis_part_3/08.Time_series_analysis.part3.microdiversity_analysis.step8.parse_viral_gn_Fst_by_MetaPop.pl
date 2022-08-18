#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum);

# Aim: Parse the viral species Fst from "MetaPop.for_each_year" folder

# Step 1 Store the Fst result
my %Viral_scf2year_pair2fst = (); # $viral_scf => $year_pair => $fst
my %Year_pair = (); # $year_pair => 1
open IN, "MetaPop.for_each_year/MetaPop/10.Microdiversity/fixation_index.tsv";
while (<IN>){
	chomp;
	if (!/^row/){
		my @tmp = split (/\t/);
		my $row_samp = $tmp[0];
		my $col_samp = $tmp[1];
		my $viral_scf = $tmp[2];
		my $fst = $tmp[3];
		
		my @Year_pair = ();
		my ($year_1) = $row_samp =~ /^(\d+?)\.viral/;
		my ($year_2) = $col_samp =~ /^(\d+?)\.viral/;
		
		push @Year_pair, $year_1;
		push @Year_pair, $year_2;
		
		@Year_pair = sort @Year_pair;
		my $year_pair = join("\_", @Year_pair);
		$Year_pair{$year_pair} = 1;
		
		if ($fst and $fst ne "NA" and $fst >= 0.01){
			$Viral_scf2year_pair2fst{$viral_scf}{$year_pair} = $fst;
		}
	}
}
close IN;

# Step 2 Get high Fst result
my %Viral_scf2year_pair2fst_high = (); # $viral_scf => $year_pair => $fst
my %Year_pair_for_fst_high = ();
foreach my $viral_scf (sort keys %Viral_scf2year_pair2fst){
	foreach my $year_pair (sort keys %Year_pair){
		if (exists $Viral_scf2year_pair2fst{$viral_scf}{$year_pair}){
			my $fst = $Viral_scf2year_pair2fst{$viral_scf}{$year_pair};
			if ($fst and $fst ne "NA" and $fst >= 0.5){
				$Viral_scf2year_pair2fst_high{$viral_scf}{$year_pair} = $fst;
				$Year_pair_for_fst_high{$year_pair} = 1;
			}
		}
	}
}

## Write down high Fst result
open OUT, ">MetaPop.for_each_year/MetaPop/Viral_scf2year_pair2fst_high.txt";
my $row=join("\t", sort keys %Year_pair_for_fst_high);
print OUT "Head\t$row\n";
foreach my $tmp1 (sort keys %Viral_scf2year_pair2fst_high){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Year_pair_for_fst_high) {       
                if (exists $Viral_scf2year_pair2fst_high{$tmp1}{$tmp2}){
                        push @tmp, $Viral_scf2year_pair2fst_high{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"NA";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 2 Get high Fst result between 2000-2003 and 2016-2019 period
my %Viral_scf2year_pair_between_2000_2003_and_2016_2019_period2fst_high = (); # $viral_scf => $year_pair => $fst_high
my %Year_pair_between_2000_2003_and_2016_2019_for_fst_high = ();
foreach my $viral_scf (sort keys %Viral_scf2year_pair2fst_high){
	foreach my $year_pair (sort keys %Year_pair_for_fst_high){
		my ($year_1, $year_2) = $year_pair =~ /^(\d+?)\_(\d+?)$/;
		if (($year_1 eq "2000" or $year_1 eq "2001" or $year_1 eq "2002" or $year_1 eq "2003") and ($year_2 eq "2016" or $year_2 eq "2017" or $year_2 eq "2018" or $year_2 eq "2019")){
			$Year_pair_between_2000_2003_and_2016_2019_for_fst_high{$year_pair} = 1;
			if (exists $Viral_scf2year_pair2fst_high{$viral_scf}{$year_pair}){
				my $fst_high = $Viral_scf2year_pair2fst_high{$viral_scf}{$year_pair};
				$Viral_scf2year_pair_between_2000_2003_and_2016_2019_period2fst_high{$viral_scf}{$year_pair} = $fst_high;
			}
		}
	}
}

## Write down high Fst result between 2000-2003 and 2016-2019 period
open OUT, ">MetaPop.for_each_year/MetaPop/Viral_scf2year_pair_between_2000_2003_and_2016_2019_period2fst_high.txt";
my $row2=join("\t", sort keys %Year_pair_between_2000_2003_and_2016_2019_for_fst_high);
print OUT "Head\t$row2\n";
foreach my $tmp1 (sort keys %Viral_scf2year_pair_between_2000_2003_and_2016_2019_period2fst_high){
        print OUT $tmp1."\t";
        my @tmp = ();
        foreach my $tmp2 (sort keys %Year_pair_between_2000_2003_and_2016_2019_for_fst_high) {       
                if (exists $Viral_scf2year_pair_between_2000_2003_and_2016_2019_period2fst_high{$tmp1}{$tmp2}){
                        push @tmp, $Viral_scf2year_pair_between_2000_2003_and_2016_2019_period2fst_high{$tmp1}{$tmp2};
                }else{              
                        push @tmp,"NA";
                }
        }
        print OUT join("\t",@tmp)."\n";
}
close OUT;

# Step 3 Get high Fst result between 2000-2003 and 2016-2019 period for viral species containing four AMGs
## Step 3.1 Store the viral species containing four AMGs
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

## Step 3.2 Write down high Fst result between 2000-2003 and 2016-2019 period for viral species containing four AMGs
open OUT, ">MetaPop.for_each_year/MetaPop/Viral_scf2year_pair_between_2000_2003_and_2016_2019_period2fst_high.for_four_AMGs.txt";
my $row3=join("\t", sort keys %Year_pair_between_2000_2003_and_2016_2019_for_fst_high);
print OUT "Head\t$row3\n";
foreach my $tmp1 (sort keys %Viral_scf2year_pair_between_2000_2003_and_2016_2019_period2fst_high){
        my ($viral_gn) = $tmp1 =~ /^(.+?)\_\_Ga/;
		if (exists $Viral_species_containing_four_AMGs{$viral_gn}){
			print OUT $tmp1."\t";
			my @tmp = ();
			foreach my $tmp2 (sort keys %Year_pair_between_2000_2003_and_2016_2019_for_fst_high) {       
					if (exists $Viral_scf2year_pair_between_2000_2003_and_2016_2019_period2fst_high{$tmp1}{$tmp2}){
							push @tmp, $Viral_scf2year_pair_between_2000_2003_and_2016_2019_period2fst_high{$tmp1}{$tmp2};
					}else{              
							push @tmp,"NA";
					}
			}
			print OUT join("\t",@tmp)."\n";
		}
}
close OUT;
