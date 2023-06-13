#!/usr/bin/perl

use strict;
use warnings;

# AIM: Find prophage hosts through scaffold-to-MAG affiliation

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
open IN, "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.txt";
while (<IN>){
	chomp;
	if (!/^tymeflies/){
		my @tmp = split (/\t/);
		my $mag = $tmp[5];
		my $num_in_cluster = $tmp[15];		
		if ($num_in_cluster ne "NA"){
			my ($img) = $mag =~ /_(33\d+?)_/;
			my @Contigs = (); # Store all the contigs into an array
			my $MAG_addr = "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/".$img."/".$mag.".fasta";
			my %MAG_seq = _store_seq("$MAG_addr");
			foreach my $header (sort keys %MAG_seq){
				my ($contig) = $header =~ /^>(.+?)$/;
				push @Contigs, $contig;	
			}			
			my $lineage = join(";", @tmp[16..22]);
			my $scaffolds = join(',', @Contigs); # Store all scaffolds in each MAG
			$TYMEFLIES_MAG_stat{$mag}[0] = $lineage;
			$TYMEFLIES_MAG_stat{$mag}[1] = $scaffolds;
		}
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

# Step 4. Write down result
open OUT, ">Host_prediction/Prophage_gn2host.txt";
foreach my $prophage_gn (sort keys %Prophage_gn2host){
	print OUT "$prophage_gn\t$Prophage_gn2host{$prophage_gn}[0]\t$Prophage_gn2host{$prophage_gn}[1]\n";
}
close OUT;



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