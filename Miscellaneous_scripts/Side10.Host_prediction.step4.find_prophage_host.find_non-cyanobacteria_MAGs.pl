#!/usr/bin/perl

use strict;
use warnings;

# AIM: Get Viral_gn_with_psbAD2host_mag.by_prophage_host result
#      viral_gn_with_psbAD => host_mags by the method of prophage host


# Step 1. Make scf to phage genome map
my %Scf2phage_gn = ();
open IN, "find /storage1/data11/TYMEFLIES_phage/*/VIBRANT_*.a.v2.min2000 -name 'phage_results.min2000.txt' |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($img_id) = $file =~ /\/storage1\/data11\/TYMEFLIES_phage\/(.+?)\//;
	my %Prophage_scf = (); # $prophage_scf => 1
	                       # Store all the prophage scaffolds in this metagenome; 
						   # for example: Ga0453107_0000053_fragment_1 => 1
	open INN, "$file";
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		if ($tmp[1] eq "prophage"){
			$Prophage_scf{$tmp[0]} = 1;
		}
	}
	close INN;
	
	`grep 'fragment' /storage1/data11/TYMEFLIES_phage/$img_id/vRhyme_best_bins_fasta_parsed/*.fasta > /storage1/data11/TYMEFLIES_phage/$img_id/vRhyme_best_bins_fasta_parsed/tmp.fragment.header.txt`;
	
	open INN, "/storage1/data11/TYMEFLIES_phage/$img_id/vRhyme_best_bins_fasta_parsed/tmp.fragment.header.txt";
	while (<INN>){
		chomp;
		my $line = $_;
		my ($prophage_scf) = $line =~ /vRhyme_unbinned\d+?__(.+?)$/;
		if (exists $Prophage_scf{$prophage_scf}){
			my ($prophage_gn) = $line =~ /\/vRhyme_best_bins_fasta_parsed\/(.+?)\.fasta\:/;
			my ($scf) = $prophage_scf =~ /^(.+?)\_fragment/;
			$Scf2phage_gn{$scf} = $prophage_gn;
		}
	}
	close INN;
	
	`rm /storage1/data11/TYMEFLIES_phage/$img_id/vRhyme_best_bins_fasta_parsed/tmp.fragment.header.txt`;
}
close IN;

# Step 2. Store TYMEFLIES MAG stat
my %TYMEFLIES_MAG_stat = (); # $mag => [0] GTDB tax [1] scaffolds
open IN, "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.CheckM_passed.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $scfs = $tmp[4];
		my $tax = $tmp[1];
		$TYMEFLIES_MAG_stat{$mag}[0] = $tax;
		$TYMEFLIES_MAG_stat{$mag}[1] = $scfs;
	}
}
close IN;

# Step 3. Find prophage hosts
my %Prophage_gn2host = (); # $prophage_gn => [0] $mag [1] $gtdb_tax
foreach my $mag (sort keys %TYMEFLIES_MAG_stat){
	my $gtdb_tax = $TYMEFLIES_MAG_stat{$mag}[0];
	my $scaffolds = $TYMEFLIES_MAG_stat{$mag}[1];
	my @Scaffolds = split (/\,/,$scaffolds);
	
	my $prophage_gn = "";
	foreach my $scf (@Scaffolds){
		if (exists $Scf2phage_gn{$scf}){
			$prophage_gn = $Scf2phage_gn{$scf};
			$Prophage_gn2host{$prophage_gn}[0] = $mag;
			$Prophage_gn2host{$prophage_gn}[1] = $gtdb_tax;
		}
	}
}

# Step 4 Find non-cyanobacteria MAG hosts that are corresponding to 
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
	my %MAGs_all = (); # Store all the MAGs 
	if (exists $Prophage_gn2host{$viral_gn}[0]){
		my $mag = $Prophage_gn2host{$viral_gn}[0];
		$MAGs_all{$mag} = 1;
	}
	
	if (%MAGs_all){
		my @MAGs = sort keys %MAGs_all;
		my $mags = join("\,", @MAGs);
		my @Taxs = ();
		foreach my $mag (@MAGs){
			my $tax = $TYMEFLIES_MAG_stat{$mag}[0];
			push @Taxs, $tax;
		}
		my $taxs = join("\,", @Taxs);
		
		$Viral_gn_with_psbAD2host_mag{$viral_gn} = $mags."\|".$taxs;
	}	
}

## Write down the result
open OUT, ">Viral_gn_with_psbAD2host_mag.by_prophage_host.txt";
foreach my $viral_gn (sort keys %Viral_gn_with_psbAD2host_mag){
	print OUT "$viral_gn\t$Viral_gn_with_psbAD2host_mag{$viral_gn}\n";
}
close OUT;
