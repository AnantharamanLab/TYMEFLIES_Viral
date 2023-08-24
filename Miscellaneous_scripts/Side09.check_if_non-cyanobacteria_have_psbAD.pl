#!/usr/bin/env perl

use strict;
use warnings;

# Aim: Check if non-cyanobacteria have psbAD
# Note: psbA: K02703 and psbD: K02706

# Step 1 Store TYMEFLIES bin summary file
my %MAG2tax = (); # $mag => $tax
my %Scf2MAG = (); # $scf => $mag
open IN, "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.CheckM_passed.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $tax = $tmp[1];
		my $scfs = $tmp[4];
		$MAG2tax{$mag} = $tax;
		my @Scfs = split("\,", $scfs);
		foreach my $scf (@Scfs){
			$Scf2MAG{$scf} = $mag
		}
	}
}
close IN;

# Step 2 Store all MAG KO hits
## Step 2.1 Store all contig to MAG info and MAG to tax info
my %TYMEFLIES_contig2MAG = (); # $contig => $mag
open IN, "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.CheckM_passed.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $contigs = $tmp[4];
		my $tax = $tmp[1];
		my @Contigs = split (/\,/, $contigs);
		foreach my $contig (@Contigs){
			$TYMEFLIES_contig2MAG{$contig} = $mag;
		}
	}
}
close IN;

my %Contig2MAG = %TYMEFLIES_contig2MAG;

## Step 2.2 Store all contigs that are determined to be viral sequence mis-binning
my %Contig_viral_seq_in_MAGs = (); # $contig => $mag (the MAG that contains the contig)
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
	}
	close INN;
}
close IN;

### Step 2.3 Store KO2cutoff_and_type hash
my %KO2cutoff_and_type = (); # $ko => $cutoff_and_type
open IN, "/slowdata/data1/Genome_profile_software/kofam_database/ko_list";
while (<IN>){
	chomp;
	if (!/^knum/){
		my @tmp = split (/\t/);
		my $ko = $tmp[0];
		my $cutoff_and_type = $tmp[1]."\t".$tmp[2];
		$KO2cutoff_and_type{$ko} = $cutoff_and_type;
	}
}
close IN;

## Step 2.4 Read hmmsearch result
my %KO2microbial_pros = (); # $ko => $microbial_pros (collection of $microbial_pro, separated by "\,")
my %MAG2kos = (); # $mag => $kos (collection of $ko, separated by "\t")
open IN, "find Host_prediction/Find_host_based_on_AMG/ -name '*.hmmsearch_result.txt' | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($ko) = $file =~ /^.+\/(K.+?)\.hmmsearch_result\.txt/;
	my ($cutoff, $type) = $KO2cutoff_and_type{$ko} =~ /^(.+?)\t(.+?)$/;

	if ($cutoff eq "\-"){
		$cutoff = 50;
	}	
	
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^#/){
			my $line = $_;
			$line =~ s/ +/ /g;
			my @tmp = split (/\s/, $line);
			my $microbial_pro = $tmp[0];
			my $score_full = $tmp[5];
			my $score_domain = $tmp[8];
			my ($microbial_scf) = $microbial_pro =~ /^(.+)\_\d+?$/;
			# Exclude all contigs/scaffolds that are determined to be viral sequence mis-binning
			if (! exists $Contig_viral_seq_in_MAGs{$microbial_scf}){ 
				if ($type eq "full"){
					if ($score_full >= $cutoff){						
						my $mag = $Scf2MAG{$microbial_scf};
						if (! exists $MAG2kos{$mag}){
							$MAG2kos{$mag} = $ko;
						}else{
							$MAG2kos{$mag} .= "\t".$ko;
						}
					}
				}else{
					if ($score_domain >= $cutoff){
						my $mag = $Scf2MAG{$microbial_scf};
						if (! exists $MAG2kos{$mag}){
							$MAG2kos{$mag} = $ko;
						}else{
							$MAG2kos{$mag} .= "\t".$ko;
						}
					}
				}
			}
		}
	}
	close INN;
}
close IN;

# Step 3 Make the %PsbA_mag_n_tax and %PsbD_mag_n_tax
my %PsbA_mag_n_tax = (); # $mag => $tax ($mag contains psbA K02703 hits)
my %PsbD_mag_n_tax = (); # $mag => $tax ($mag contains psbD K02706 hits)
foreach my $mag (sort keys %MAG2kos){
	my $kos = $MAG2kos{$mag};
	if ($kos =~ /K02703/){
		$PsbA_mag_n_tax{$mag} = $MAG2tax{$mag};
	}elsif ($kos =~ /K02706/){
		$PsbD_mag_n_tax{$mag} = $MAG2tax{$mag};
	}	
}

# Step 4 Write down both hash
open OUT, ">PsbA_mag_n_tax.txt";
foreach my $key (sort keys %PsbA_mag_n_tax){
	print OUT "$key\t$PsbA_mag_n_tax{$key}\n";
}
close OUT;

open OUT, ">PsbD_mag_n_tax.txt";
foreach my $key (sort keys %PsbD_mag_n_tax){
	print OUT "$key\t$PsbD_mag_n_tax{$key}\n";
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