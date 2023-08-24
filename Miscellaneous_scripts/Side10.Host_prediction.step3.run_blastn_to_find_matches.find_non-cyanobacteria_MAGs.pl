#!/usr/bin/perl

use strict;
use warnings;

# AIM: Get Viral_gn_with_psbAD2host_mag.by_crispr_match result
#      viral_gn_with_psbAD => host_mags by the method of crispr match searching


# Step 1. Parse blastn results to store matches
# There are two categories of matches: 
# (1) have 0 or 1 mismatch over the entire spacer length (‘CRISPR (near)identical’) 
# (2) have ≥80% identity over the entire spacer length (‘CRISPR multiple partial’)

## Step 1.1 Store the spacer length
my %Spacer2length = (); # $spacer => $length; for example: 2004178001.a:gws2_d1_0103_30__Spacer_0001 => 32 
my %Hash = _store_seq("Host_prediction/All_spacers.fasta");
foreach my $key (sort keys %Hash){
	my ($key_clean) = $key =~ /^>(.+?)$/;
	my $length = length($Hash{$key});
	$Spacer2length{$key_clean} = $length;
}

## Step 1.2 Parse blastn results
my %CRISPR_near_identical_matches = (); # $match => 1; an example for a match: 2004178001.a:gws2_d1_0103_30__Spacer_0006|3300044846__vRhyme_449__Ga0453141_0000520
my %CRISPR_partial_matches = (); # $match => 1;
open IN, "Host_prediction/all_phage_genome2all_spacer.blastn_out.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $phage_scf = $tmp[0];
	my $spacer = $tmp[1];
	
	my $spacer_length = $Spacer2length{$spacer};
	
	my $identity = $tmp[2];
	my $mismatches = $tmp[4];
	my $gap_openings = $tmp[5];
	
	my $s_start = $tmp[8];
	my $s_end = $tmp[9];
	my $spacer_matched_length = abs($s_end - $s_start) + 1;
	
	if ($spacer_matched_length eq $spacer_length){
		if ($mismatches <= 1 and $gap_openings == 0){
			$CRISPR_near_identical_matches{"$phage_scf\|$spacer"} = 1;
		}elsif ($identity >= 90 and !($mismatches <= 1 and $gap_openings == 0)){
			$CRISPR_partial_matches{"$phage_scf\|$spacer"} = 1;
		}
	}
}
close IN;

my $num_CRISPR_near_identical_matches = scalar (keys (%CRISPR_near_identical_matches));
print "There are $num_CRISPR_near_identical_matches CRISPR (near)identical matches\n";
my $num_CRISPR_partial_matches = scalar (keys (%CRISPR_partial_matches));
print "There are $num_CRISPR_partial_matches CRISPR partial matches\n";

# Step 2. Predict host taxonomy
## Step 2.1 Store phage scf to gn hash
my %Phage_scf2gn = (); # $scf => $gn
my %Phage_gn2scf = (); # $gn => $scf collection separated by "\t"
open IN, "Host_prediction/All_phage_genomes.fasta";
while (<IN>){
	chomp;
	if (/^>/){
		my $line = $_;
		my ($scf) = $line =~ /^>(.+?)$/;
		my ($gn) = $scf =~ /^(.+?\_\_.+?)\_\_.+?$/;
		$Phage_scf2gn{$scf} = $gn;
		if (!exists $Phage_gn2scf{$gn}){
			$Phage_gn2scf{$gn} = $scf;
		}else{
			$Phage_gn2scf{$gn} .= "\t".$scf;
		}
	}
}
close IN;

## Step 2.2 Store spacer to taxonomy hash
my %TYMEFLIES_MAG2scfs = (); # $mag => $scf joined by ','
my %MAG_scf2lineage = (); # $scf => $lineage
open IN, "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.CheckM_passed.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $scfs = $tmp[4];
		my $lineage = $tmp[1];
		$TYMEFLIES_MAG2scfs{$mag} = $scfs;
		
		my @Scfs = split (/\,/, $scfs);
		foreach my $scf (@Scfs){
			$MAG_scf2lineage{$scf} = $lineage;
		}
	
	}
}
close IN;

my %Spacer2lineage = (); # $spacer => $lineage
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/All_spacers.fasta";
while (<IN>){
	chomp;
	if (/^>/){
		my $line = $_;
		my ($spacer) = $line =~ /^>(.+?)$/;
		my ($scf) = $spacer =~ /^(.+?)\_\_Spacer/;
		my $lineage = "NA";
		if (exists $MAG_scf2lineage{$scf}){
			$lineage = $MAG_scf2lineage{$scf};	
		}
		$Spacer2lineage{$spacer} = $lineage;
	}
}
close IN;

## Print %Spacer2lineage
open OUT, ">Host_prediction/Spacer2lineage.txt";
foreach my $spacer (sort keys %Spacer2lineage){
	print OUT "$spacer\t$Spacer2lineage{$spacer}\n";
}
close OUT;

## Step 2.3 Store phage gn to all lineages (corresponding to its scaffolds)
my %Phage_gn2lineage = (); # $gn => [0] lineages (separated by "\t") based on CRISPR_near_identical_matches
                                   #[1] lineages (separated by "\t") based on CRISPR_partial_matches
my %Phage_gn2spacer_based_on_CRISPR_near_identical_matches = (); # $gn => $spacers (collection of $spacer, separated by "\,")
my %Phage_gn2spacer_based_on_CRISPR_partial_matches = (); # $gn => $spacers (collection of $spacer, separated by "\,")
foreach my $match (sort keys %CRISPR_near_identical_matches){
	my ($phage_scf,$spacer) = $match =~ /^(.+?)\|(.+?)$/;
	my $phage_gn = $Phage_scf2gn{$phage_scf};
	my $lineage = "";
	if (exists $Spacer2lineage{$spacer}){
		$lineage = $Spacer2lineage{$spacer};
	}
	
	if ($lineage and $lineage ne 'NA'){ # Do not store "NA" lineage
		if (!exists $Phage_gn2lineage{$phage_gn}[0]){
			$Phage_gn2lineage{$phage_gn}[0] = $lineage;
		}else{
			$Phage_gn2lineage{$phage_gn}[0] .= "\t".$lineage;
		}
	}
	
	if (!exists $Phage_gn2spacer_based_on_CRISPR_near_identical_matches{$phage_gn}){
		$Phage_gn2spacer_based_on_CRISPR_near_identical_matches{$phage_gn} = $spacer;
	}else{
		$Phage_gn2spacer_based_on_CRISPR_near_identical_matches{$phage_gn} .= "\,".$spacer;
	}
}	

foreach my $match (sort keys %CRISPR_partial_matches){
	my ($phage_scf,$spacer) = $match =~ /^(.+?)\|(.+?)$/;
	my $phage_gn = $Phage_scf2gn{$phage_scf};
	my $lineage = "";
	if (exists $Spacer2lineage{$spacer}){
		$lineage = $Spacer2lineage{$spacer};
	}
	
	if ($lineage and $lineage ne 'NA'){	# Do not store "NA" lineage
		if (!exists $Phage_gn2lineage{$phage_gn}[1]){
			$Phage_gn2lineage{$phage_gn}[1] = $lineage;
		}else{
			$Phage_gn2lineage{$phage_gn}[1] .= "\t".$lineage;
		}
	}
	
	if (!exists $Phage_gn2spacer_based_on_CRISPR_partial_matches{$phage_gn}){
		$Phage_gn2spacer_based_on_CRISPR_partial_matches{$phage_gn} = $spacer;
	}else{
		$Phage_gn2spacer_based_on_CRISPR_partial_matches{$phage_gn} .= "\,".$spacer;
	}
}

# Step 3 Find non-cyanobacteria MAG hosts that are corresponding to 
#        viral genome with psbAD 

## Store Viral_gn_with_psbAD_connected_to_noncyanobacteria_host.txt
my %Viral_gn_with_psbAD_connected_to_noncyanobacteria_host = ();
open IN, "Viral_gn_with_psbAD_connected_to_noncyanobacteria_host.txt";
while (<IN>){
	chomp;
	my $viral_gn = $_;
	$Viral_gn_with_psbAD_connected_to_noncyanobacteria_host{$viral_gn} = 1;
}
close IN;

my %Scf2mag = (); # $scf => $mag 
my %MAG2tax = (); # $mag => $tax 
open IN, "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.CheckM_passed.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $tax = $tmp[1];
		my $scfs = $tmp[4];
		my @Scfs = split (/\,/,$scfs);
		foreach my $scf (@Scfs){
			$Scf2mag{$scf} = $mag;
		}
		
		$MAG2tax{$mag} = $tax;
	}
}

my %Spacer2mag = (); # $spacer => $mag
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/All_spacers.fasta";
while (<IN>){
	chomp;
	if (/^>/){
		my $line = $_;
		my ($spacer) = $line =~ /^>(.+?)$/;
		my ($scf) = $spacer =~ /^(.+?)\_\_Spacer/;
		if (exists $Scf2mag{$scf}){ # Only use spacers from TYMEFLIES MAGs
			my $mag = $Scf2mag{$scf};
			$Spacer2mag{$spacer} = $mag;
		}
	}
}
close IN;

my %Viral_gn_with_psbAD2host_mag = (); # $viral_gn (with psbAD inside) => $mags (collection of $mag separated by "\,", the host) and $taxs (collection of $tax separated by "\,") 
foreach my $viral_gn (sort keys %Viral_gn_with_psbAD_connected_to_noncyanobacteria_host){
	my %MAGs_all = (); # Store all the MAGs 
	if (exists $Phage_gn2spacer_based_on_CRISPR_near_identical_matches{$viral_gn}){
		my @Spacers = split (/\,/, $Phage_gn2spacer_based_on_CRISPR_near_identical_matches{$viral_gn});
		foreach my $spacer (@Spacers){
			if (exists $Spacer2mag{$spacer}){
				my $mag = $Spacer2mag{$spacer};
				$MAGs_all{$mag} = 1;
			}
		}
	}elsif(exists $Phage_gn2spacer_based_on_CRISPR_partial_matches{$viral_gn}){
		my @Spacers = split (/\,/, $Phage_gn2spacer_based_on_CRISPR_partial_matches{$viral_gn});
		if ((scalar @Spacers) >= 2){
			foreach my $spacer (@Spacers){
				if (exists $Spacer2mag{$spacer}){
					my $mag = $Spacer2mag{$spacer};
					$MAGs_all{$mag} = 1;
				}
			}			
		}
	}
	
	if (%MAGs_all){
		my @MAGs = sort keys %MAGs_all;
		my $mags = join("\,", @MAGs);
		my @Taxs = ();
		foreach my $mag (@MAGs){
			my $tax = $MAG2tax{$mag};
			push @Taxs, $tax;
		}
		my $taxs = join("\,", @Taxs);
		
		$Viral_gn_with_psbAD2host_mag{$viral_gn} = $mags."\|".$taxs;
	}	
}

## Write down the result
open OUT, ">Viral_gn_with_psbAD2host_mag.by_crispr_match.txt";
foreach my $viral_gn (sort keys %Viral_gn_with_psbAD2host_mag){
	print OUT "$viral_gn\t$Viral_gn_with_psbAD2host_mag{$viral_gn}\n";
}
close OUT;



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
