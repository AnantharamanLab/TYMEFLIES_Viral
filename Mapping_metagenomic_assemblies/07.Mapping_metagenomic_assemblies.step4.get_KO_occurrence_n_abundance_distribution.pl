#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get the KO occurrence and abundance distribution pattern

# Step 1 Store AMG summary table
my %AMG_summary = (); # $pro => $ko
my %KOs= (); # $ko => 1;
my %IMG2date = (); # $img_id => $date_n_season
open IN, "AMG_analysis/AMG_summary.txt";
while (<IN>){
	chomp;
	if (!/^Pro/){
		my @tmp = split (/\t/);
		my $pro = $tmp[0];
		my $ko = $tmp[2];
		my $ko_detail = $tmp[3];
		my $date_n_season = $tmp[1];
		$AMG_summary{$pro} = $ko;
		my ($img_id) = $pro =~ /^(33.+?)\_/;
		$IMG2date{$img_id} = $date_n_season;
		$KOs{$ko} = $ko_detail;
	}
}
close IN;

# Step 2 Get AMG KO abundances (the metagenome size was normalized by 100M reads) in each metagenome
## Step 2.1 Store scaffold coverage for each phage genome
my %Scf2cov = (); # $scf => $cov (within individual metagenome)
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

## Step 2.2 Get AMG abundances (normalzied by 100M reads per metagenome)
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

my %KO2metagenome2abun_normalized = (); # $ko => $img_id => $abun_normalized_sum (normalized abundance sum of all AMG)
foreach my $ko (sort keys %KOs){
	foreach my $pro (sort keys %AMG_summary){
		if ($AMG_summary{$pro} eq $ko){
			my $abun_normalized = 0;
			my ($img_id,$scf) = $pro =~ /^(33\d+?)\_\_vRhyme\_.+?\_\_(Ga.+)\_\d+?$/;
			my $read_num = $IMG_ID2read_num{$img_id}; 
			my $abun = $Scf2cov{$scf}; 
			if ($abun){
				$abun_normalized = $abun / ($read_num / 100000000);
			}
			
			$KO2metagenome2abun_normalized{$ko}{$img_id} += $abun_normalized;

		}
	}
}

# Step 3 Get AMG KO occurence
my %KO2metagenomes = (); # $ko => $metagenomes (collection of $metagenome, separated by "\t")
foreach my $ko (sort keys %KOs){
	foreach my $pro (sort keys %AMG_summary){
		if ($AMG_summary{$pro} eq $ko){
			my ($metagenome) = $pro =~ /^(.+?)\_\_/;
			if (!exists $KO2metagenomes{$ko}){
				$KO2metagenomes{$ko} = $metagenome;
			}else{
				if ($KO2metagenomes{$ko} !~ /$metagenome/){
					$KO2metagenomes{$ko} .= "\t".$metagenome;
				}
			}
		
		}
	}
}

# Step 4 Store and write down %KO2occurrence_n_abundance hash 
my %KO2occurrence_n_abundance = (); # $ko => [0] $occurrence [1] $abundance
foreach my $ko (sort keys %KOs){
	my $occurrence = 0;
	my @IMG = split (/\t/, $KO2metagenomes{$ko});
	$occurrence = scalar @IMG;
	
	my $abundance = 0; # Store the mean abundance in IMG that species occur
	my @Abundance = (); # Store abundance in IMG that species occur
	foreach my $img_id (@IMG){
		my $abun = $KO2metagenome2abun_normalized{$ko}{$img_id};
		push @Abundance, $abun;
	}
	$abundance = _avg(@Abundance);
	
	$KO2occurrence_n_abundance{$ko}[0] = $occurrence;
	$KO2occurrence_n_abundance{$ko}[1] = $abundance;
}

open OUT, ">KO2occurrence_n_abundance.txt";
foreach my $ko (sort keys %KO2occurrence_n_abundance){
	print OUT "$ko\t$KO2occurrence_n_abundance{$ko}[0]\t$KO2occurrence_n_abundance{$ko}[1]\n";
}
close OUT;



# Subroutine 
sub _avg {
    my $total;
    $total += $_ foreach @_;
    # sum divided by number of components.
    return $total / @_;
}

