#!/usr/bin/perl

use strict;
use warnings;

# AIM: Get Viral_gn_with_psbAD2host_mag.by_sequence_identity result
#      viral_gn_with_psbAD => host_mags by the method of sequence searching  


# Step 1. Parse result to get viral genome to final host prediction
## Step 1.1 Store all contig to MAG info and MAG to tax info
my %TYMEFLIES_contig2MAG = (); # $contig => $mag
my %TYMEFLIES_MAG2contigs = (); # $mag => $contig collection separated by "\,"
my %MAG2tax = (); # $mag => $tax
open IN, "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.CheckM_passed.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $contigs = $tmp[4];
		my $tax = $tmp[1];
		$TYMEFLIES_MAG2contigs{$mag} = $contigs;
		
		my @Contigs = split (/\,/, $contigs);
		foreach my $contig (@Contigs){
			$TYMEFLIES_contig2MAG{$contig} = $mag;
		}
		
		$MAG2tax{$mag} = $tax;
	}
}
close IN;

my %GEM_contig2MAG = (); # $contig => $mag
my %GEM_MAG2contigs = (); # $mag => $contig collection separated by "\,"
open IN, "/storage1/databases/GEM/GEM_all_MAGs_stat.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $contigs = $tmp[4];
		my $tax = $tmp[1];
		$GEM_MAG2contigs{$mag} = $contigs;
		
		my @Contigs = split (/\,/, $contigs);
		foreach my $contig (@Contigs){
			$GEM_contig2MAG{$contig} = $mag;
		}
		
		$MAG2tax{$mag} = $tax;		
	}
}
close IN;

my %Contig2MAG = (%TYMEFLIES_contig2MAG, %GEM_contig2MAG);
my %MAG2contigs = (%TYMEFLIES_MAG2contigs, %GEM_MAG2contigs);

## Step 1.2 Store all contigs that are determined to be viral sequence mis-binning,
#           and store the match result of viral sequence to contig   
my %Contig_viral_seq_in_MAGs = (); # $contig => $mag (the MAG that contains the contig)
my %Viral_seq2contig = (); # $viral_seq => $contig collection separated by "\,"
open IN, "ls /storage1/data11/TYMEFLIES_phage/Robin_MAGs/split_fsa/*.fsa |";
while (<IN>){
	chomp;
	my $fsa_file = $_;
	my ($out_name) = $fsa_file =~ /split_fsa\/(.+?)\.fsa/;
	my %Seq = _store_seq("$fsa_file");
	my %Contig2len = (); # $contig => $len
	foreach my $key (sort keys %Seq){
		my ($key_clean) = $key =~ /^>(.+?)$/;
		my $len = length($Seq{$key});
		$Contig2len{$key_clean} = $len;
	}
	
	# Read the blastn out
	open INN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/filter_MAGs_out/TYMEFLIES_MAGs_${out_name}.blastn_out.txt";
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $contig = $tmp[0];
		my $viral_seq = $tmp[1];
		my $qstart = $tmp[6];
		my $qend = $tmp[7];
		my $tstart = $tmp[8];		
		my $tend = $tmp[9];		
		my $iden = $tmp[2];
		my $host_covered_len = abs($qend - $qstart);
		my $viral_covered_len = abs($tend - $tstart);		
		my $host_covered_prec = $host_covered_len / $Contig2len{$contig};
		if ($host_covered_prec >= 0.5 and $iden >= 80){ # Cutoff: viral sequence hit covered region ≥ 50% of the host contig
			my $mag = $Contig2MAG{$contig};
			$Contig_viral_seq_in_MAGs{$contig} = $mag;
		}
		if ($host_covered_len >= 2000 and $viral_covered_len >= 2000 and $iden >= 95){ #  Host predictions were then based on matches of ≥90% nucleotide identity covering ≥2 kb of the virus and (putative) host sequences.
			if (! exists $Viral_seq2contig{$viral_seq}){
				$Viral_seq2contig{$viral_seq} = $contig;
			}else{
				$Viral_seq2contig{$viral_seq} .= "\,".$contig;
			}
		}
	}
	close INN;
}
close IN;

open IN, "ls /storage1/databases/GEM/split_fsa/*.fsa |";
while (<IN>){
	chomp;
	my $fsa_file = $_;
	my ($out_name) = $fsa_file =~ /split_fsa\/(.+?)\.fsa/;
	my %Seq = _store_seq("$fsa_file");
	my %Contig2len = (); # $contig => $len
	foreach my $key (sort keys %Seq){
		my ($key_clean) = $key =~ /^>(.+?)$/;
		my $len = length($Seq{$key});
		$Contig2len{$key_clean} = $len;
	}
	
	# Read the blastn out
	open INN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/filter_MAGs_out/GEM_MAGs_${out_name}.blastn_out.txt";
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $contig = $tmp[0];
		my $viral_seq = $tmp[1];
		my $qstart = $tmp[6];
		my $qend = $tmp[7];
		my $tstart = $tmp[8];		
		my $tend = $tmp[9];		
		my $iden = $tmp[2];
		my $host_covered_len = abs($qend - $qstart);
		my $viral_covered_len = abs($tend - $tstart);		
		my $host_covered_prec = $host_covered_len / $Contig2len{$contig};
		if ($host_covered_prec >= 0.5 and $iden >= 80){ # Cutoff: viral sequence hit covered region ≥ 50% of the host contig
			my $mag = $Contig2MAG{$contig};
			$Contig_viral_seq_in_MAGs{$contig} = $mag;
		}
		if ($host_covered_len >= 2000 and $viral_covered_len >= 2000 and $iden >= 95){ #  Host predictions were then based on matches of ≥90% nucleotide identity covering ≥2 kb of the virus and (putative) host sequences.
			if (! exists $Viral_seq2contig{$viral_seq}){
				$Viral_seq2contig{$viral_seq} = $contig;
			}else{
				$Viral_seq2contig{$viral_seq} .= "\,".$contig;
			}
		}
	}
	close INN;
}
close IN;

# Step 2 Find non-cyanobacteria MAG hosts that are corresponding to 
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

my %Viral_gn_with_psbAD2host_mag = (); # $viral_gn (with psbAD inside) => $mags (collection of $mag separated by "\,", the host) and $taxs (collection of $tax separated by "\,") 
foreach my $viral_gn (sort keys %Viral_gn_with_psbAD_connected_to_noncyanobacteria_host){
	my %Contigs_all = (); # Store the corresponding MAG contigs
	foreach my $viral_seq (sort keys %Viral_seq2contig){
		if ($viral_seq =~ /__vRhyme_/){ # Only use the viral seq from TYMEFLIES project (contains "__vRhyme_" inside)
			my ($viral_gn2) = $viral_seq =~ /^(.+?\_\_.+?)\_\_/;
			if ($viral_gn2 eq $viral_gn){
				my @Contigs = split (/\,/, $Viral_seq2contig{$viral_seq});
				foreach my $contig (@Contigs){
					if (! exists $Contig_viral_seq_in_MAGs{$contig}){ # To exclude contigs that are potentially viral sequences
						$Contigs_all{$contig} = 1;
					}
				}
			}
		}
	}
	
	my %MAGs_all = (); # Store the corresponding MAGs 
	if (%Contigs_all){
		foreach my $contig (sort keys %Contigs_all){
			my $mag = $Contig2MAG{$contig};
			if (exists $TYMEFLIES_MAG2contigs{$mag}){
				$MAGs_all{$mag} = 1; # Only use the genome from TYMEFLIES
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
open OUT, ">Viral_gn_with_psbAD2host_mag.by_sequence_identity.txt";
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
