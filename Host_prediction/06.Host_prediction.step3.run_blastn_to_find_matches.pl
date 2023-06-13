#!/usr/bin/perl

use strict;
use warnings;

# AIM: Find MAG crispr spacer to phage matches by blastn

=pod
# Step 1. Get all the spacers from TYMEFLIES MAGs
## Step 1.1 Store all MAG ID in the folder "Host_prediction/find_crispr_out"
my %MAG_ID = (); # Store all the MAG ID; $mag_id => 1; for example: TYMEFLIES_MAGs_ME2015-07-03_3300042555_group6_bin60 => 1
open IN, "find /storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out -name '*.minced_out.crisprs' |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($mag_id) = $file =~ /find_crispr_out\/(.+?)\.minced_out\.crisprs/;
	$MAG_ID{$mag_id} = 1;
}
close IN;

## Step 1.2 Store %Seq2spacers for MAG in a loop
my %Seq2spacers = (); # Store all spacers to seqs; # $seq => $spacers (combining both MinCED and PILER-CR results)
foreach my $img_id (sort keys %MAG_ID){
	my %Seq2spacers_this_MAG = (); # $seq => $spacers
	my %Seq2spacers_MinCED = ();
	if (-s "/storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.minced_out.crisprs"){
		%Seq2spacers_MinCED = _store_minced_crispr_out("/storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.minced_out.crisprs");
	}
	
	my %Seq2spacers_PILER_CR = ();
	`sed -n '/DETAIL REPORT/,/SUMMARY BY SIMILARITY/p' /storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.pilercr_out.crisprs > /storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.pilercr_out.crisprs.parsed`;
	if (-s "/storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.pilercr_out.crisprs.parsed"){
		%Seq2spacers_PILER_CR = _store_pilercr_crispr_out("/storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.pilercr_out.crisprs.parsed");
	}
	`rm /storage1/data11/TYMEFLIES_phage/Host_prediction/find_crispr_out/$img_id.pilercr_out.crisprs.parsed`;
	
	if (!(%Seq2spacers_MinCED) and !(%Seq2spacers_PILER_CR)){
		%Seq2spacers_this_MAG = ();
	}elsif (!(%Seq2spacers_MinCED) and %Seq2spacers_PILER_CR){
		%Seq2spacers_this_MAG =  %Seq2spacers_PILER_CR;
	}elsif (%Seq2spacers_MinCED and !(%Seq2spacers_PILER_CR)){
		%Seq2spacers_this_MAG = %Seq2spacers_MinCED;
	}elsif (%Seq2spacers_MinCED and %Seq2spacers_PILER_CR){
		%Seq2spacers_this_MAG = %Seq2spacers_MinCED;
		foreach my $seq (sort keys %Seq2spacers_this_MAG){
			if (exists $Seq2spacers_PILER_CR{$seq}){
				my $spacers_MinCED = $Seq2spacers_MinCED{$seq};
				my $spacers_PILER_CR = $Seq2spacers_PILER_CR{$seq};
				my $spacers_combined = $spacers_MinCED."\,".$spacers_PILER_CR;
				my @Spacers_combined = split (/\,/,$spacers_combined);
				my %Spacers_combined_drep = ();
				foreach my $key (@Spacers_combined){
					$Spacers_combined_drep{$key} = 1;
				}				
				$spacers_combined = join("\,", (keys %Spacers_combined_drep)); # dRep-ed spacers
				$Seq2spacers_this_MAG{$seq} = $spacers_combined;
			}
		}
		
		foreach my $seq (sort keys %Seq2spacers_PILER_CR){
			if (!exists $Seq2spacers_this_MAG{$seq}){
				$Seq2spacers_this_MAG{$seq} = $Seq2spacers_PILER_CR{$seq};
			}
		}
	}
	
	%Seq2spacers = (%Seq2spacers, %Seq2spacers_this_MAG);
}

## Step 1.3 Write down all the spacers into a fasta files
my %All_spacers = (); # $head => $spacer_seq
foreach my $seq (sort keys %Seq2spacers){
	my $spacers = $Seq2spacers{$seq};
	my @Spacers = split (/\,/,$spacers);
	for(my $i=0; $i<=$#Spacers; $i++){
		my $j = $i+1;
		$j = (sprintf "%04d", $j);
		my $head = ">$seq\_\_Spacer_$j";
		$All_spacers{$head} = $Spacers[$i];
	}
}

open OUT, ">/storage1/data11/TYMEFLIES_phage/Host_prediction/All_spacers.fasta";
foreach my $key (sort keys %All_spacers){
	print OUT "$key\n$All_spacers{$key}\n";
}
close OUT;

# Step 2. Make blastn db 

`makeblastdb -in Host_prediction/All_spacers.fasta -title All_spacers -dbtype nucl -out Host_prediction/All_spacers_blastndb`;

# Step 3. Run blastn to find spacer and phage association
`find /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed -name '*.fasta' -exec cat {} + > /storage1/data11/TYMEFLIES_phage/Host_prediction/All_phage_genomes.fasta`;

`mkdir Host_prediction/All_phage_genomes_split_fsa`;
`perl /storage1/data11/TYMEFLIES_phage/split_multifasta.pl --in Host_prediction/All_phage_genomes.fasta --output_dir=Host_prediction/All_phage_genomes_split_fsa --seqs_per_file=10000`;

open OUT, ">tmp.run_blastn_to_find_matches_in_parallel.sh";
open IN, "ls Host_prediction/All_phage_genomes_split_fsa/*.fsa |";
while (<IN>){
	chomp;
	my $fsa = $_;
	my ($fsa_name) = $fsa =~ /All_phage_genomes_split_fsa\/(.+?)\.fsa/;
	print OUT "blastn -task blastn-short -evalue 1 -dust no -word_size 7 -num_threads 1 -outfmt 6 -query $fsa -db Host_prediction/All_spacers_blastndb -out Host_prediction/all_phage_genome.${fsa_name}.2all_spacer.blastn_out.txt\n";
}
close IN;
close OUT;

`cat tmp.run_blastn_to_find_matches_in_parallel.sh | parallel -j 20`;
`rm tmp.run_blastn_to_find_matches_in_parallel.sh`;
=cut

# Step 4. Parse blastn results to store matches
# There are two categories of matches: 
# (1) have 0 or 1 mismatch over the entire spacer length (‘CRISPR (near)identical’) 
# (2) have ≥ 90% identity over the entire spacer length (‘CRISPR multiple partial’)

## Step 4.1 Store the spacer length
my %Spacer2length = (); # $spacer => $length; for example: 2004178001.a:gws2_d1_0103_30__Spacer_0001 => 32 
my %Hash = _store_seq("Host_prediction/All_spacers.fasta");
foreach my $key (sort keys %Hash){
	my ($key_clean) = $key =~ /^>(.+?)$/;
	my $length = length($Hash{$key});
	$Spacer2length{$key_clean} = $length;
}

## Step 4.2 Parse blastn results

`cat Host_prediction/all_phage_genome.*.2all_spacer.blastn_out.txt > Host_prediction/all_phage_genome2all_spacer.blastn_out.txt`;
`rm Host_prediction/all_phage_genome.*.2all_spacer.blastn_out.txt`;

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

# Step 5. Predict host taxonomy
## Step 5.1 Store phage scf to gn hash
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

## Step 5.2 Store spacer to taxonomy hash
my %TYMEFLIES_MAG2scfs = (); # $mag => $scf joined by ','
my %MAG_scf2lineage = (); # $scf => $lineage
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
			$TYMEFLIES_MAG2scfs{$mag} = join(',', @Contigs); # Store all scaffolds in each MAG
			
			foreach my $contig (@Contigs){ # Here $contig is the same with $scf
				$MAG_scf2lineage{$contig} = $lineage;
			}
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

=pod
## Get rid of TYMEFLIES MAGs that are contaminated MAGs
### Store the Non-cyanobacteria_MAGs_that_connect_with_cyanobacteria.txt (the contaminated MAGs list)
my %Non_cyanobacteria_MAGs_that_connect_with_cyanobacteria = (); # $mag => 1
open IN, "Check_non-cyanobacteria/Non-cyanobacteria_MAGs_that_connect_with_cyanobacteria.txt";
while (<IN>){
	chomp;
	my $mag = $_;
	$Non_cyanobacteria_MAGs_that_connect_with_cyanobacteria{$mag} = 1;
}
close IN;

my %TYMEFLIES_MAG_scf = (); # $scf => 1
open IN, "/storage1/data11/TYMEFLIES_phage/Binning_Data/TYMEFLIES_all_MAGs_stat.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		if (!exists $Non_cyanobacteria_MAGs_that_connect_with_cyanobacteria{$mag}){ # Get rid of those contaminated MAGs
			my $scfs = $tmp[4];
			my @Scfs = split (/\,/,$scfs);
			foreach my $scf (@Scfs){
				$TYMEFLIES_MAG_scf{$scf} = 1;
			}
		}
	}
}
close IN;

### Delete spacers that are not from TYMEFLIES MAGs
foreach my $spacer (sort keys %Spacer2lineage){
	my ($scf) = $spacer =~ /(^.+?)\_\_/;
	if (!exists $TYMEFLIES_MAG_scf{$scf}){
		delete $Spacer2lineage{$spacer}; 
	}
}
=cut
## Print %Spacer2lineage
open OUT, ">Host_prediction/Spacer2lineage.txt";
foreach my $spacer (sort keys %Spacer2lineage){
	print OUT "$spacer\t$Spacer2lineage{$spacer}\n";
}
close OUT;

## Step 5.3 Store phage gn to all lineages (corresponding to its scaffolds)
my %Phage_gn2lineage = (); # $gn => [0] lineages (separated by "\t") based on CRISPR_near_identical_matches
                                   #[1] lineages (separated by "\t") based on CRISPR_partial_matches

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
}

## Step 5.4 Store phage gn to final lineage based on 80% consensus on each rank
my %Phage_gn2lineage_final = (); # $gn => $lineage_final
foreach my $gn (sort keys %Phage_gn2lineage){
	my $lineages_by_near_identical_matches = "";
	if ($Phage_gn2lineage{$gn}[0]){
		$lineages_by_near_identical_matches = $Phage_gn2lineage{$gn}[0];
	}
	my $lineages_by_partial_matches = "";
	if ($Phage_gn2lineage{$gn}[1]){
		$lineages_by_partial_matches = $Phage_gn2lineage{$gn}[1];
	}
	
	if ($lineages_by_near_identical_matches){ # If near_identical_matches has hits
		my @tmp = split (/\t/, $lineages_by_near_identical_matches);
		my $consensus_lineage = _find_consensus_lineages_based_on_each_rank(@tmp);
		if ($consensus_lineage){
			$Phage_gn2lineage_final{$gn} = $consensus_lineage;
		}
	}elsif (! $lineages_by_near_identical_matches and $lineages_by_partial_matches){
		my @tmp = split (/\t/, $lineages_by_partial_matches); 
		if ((scalar @tmp) >= 2){ # Need to have at least two matches for "partial matches"
			my $consensus_lineage = _find_consensus_lineages_based_on_each_rank(@tmp);
			if ($consensus_lineage){
				$Phage_gn2lineage_final{$gn} = $consensus_lineage;
			}
		}
	}
}

# Step 5.5 Write down result
open OUT, ">Host_prediction/Phage_gn2host_tax_by_CRISPR_matches.txt";
foreach my $gn (sort keys %Phage_gn2lineage_final){
	print OUT "$gn\t$Phage_gn2lineage_final{$gn}\n";
}
close OUT;



# Subroutine

sub _store_minced_crispr_out{
	my $file = $_[0];
	my %Hash = (); # $spacer => $seq
	my %Seq = (); # $seq => 1
	my $seq = "";
	open IN_, $file;
	while (<IN_>){
		chomp;
		if (/^Sequence/){
			my $line = $_;
			($seq) = $line =~ /\'(.+?)\'/;
			$Seq{$seq} = 1;
		}elsif(/\[/){
			my $line = $_;
			my @tmp = split (/\t/,$line);
			my $spacer = $tmp[3];
			$Hash{$spacer} = $seq;
		}
	}
	close IN_;
	
	my %Hash2 = (); # $seq => $spacer collection separated by ","
	foreach my $spacer (sort keys %Hash){
		my $seq = $Hash{$spacer};
		if (!exists $Hash2{$seq}){
			$Hash2{$seq} = $spacer;
		}else{
			$Hash2{$seq} .= "\,".$spacer;
		}
	}
	
	return %Hash2; 
}

sub _store_pilercr_crispr_out{
	my $file = $_[0];
	my %Hash = (); # $spacer => $seq
	my %Seq = (); # $seq => 1
	my $seq = "";
	open IN_, $file;
	while (<IN_>){
		chomp;
		if (/^>/){
			my $line = $_;
			($seq) = $line =~ /^>(.+?)$/;
			$Seq{$seq} = 1;
		}elsif(/\.\./){
			my $line = $_;
			$line =~ s/ +/ /g; # Change multiple spaces into one space
			$line =~ s/^ //g; # Delete the first space in the front
			my @tmp = split (/\s/,$line);
			my $spacer_length = $tmp[3];
			if ($spacer_length =~ /\d\d/){ # Get to the spacer line
				my $spacer = $tmp[-1];
				$Hash{$spacer} = $seq;				
			}						
		}
	}
	close IN_;
	
	
	my %Hash2 = (); # $seq => $spacer collection separated by ","
	foreach my $spacer (sort keys %Hash){
		my $seq = $Hash{$spacer};
		if (!exists $Hash2{$seq}){
			$Hash2{$seq} = $spacer;
		}else{
			$Hash2{$seq} .= "\,".$spacer;
		}
	}
	
	return %Hash2; 
}

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

sub _find_consensus_lineages_based_on_each_rank{ # The default is 80% consensus
	my @Lineages = @_; # Get the passed array
	
	my %Domain = (); # $domain => the times of this domain appears
	my %Phylum = (); # $phylum => the times of this phylum appears
	my %Class = (); # $class => the times of this class appears
	my %Order = (); # $order => the times of this order appears
	my %Family = (); # $family => the times of this family appears
	my %Genus = (); # $genus => the times of this genus appears
	my %Species = (); # $species => the times of this species appears
	my $lineage_hit_num = scalar @Lineages;
	
	foreach my $lineage (@Lineages){
		my ($domain,$phylum,$class,$order,$family,$genus,$species) = $lineage =~ /^(.+?)\;(.+?)\;(.+?)\;(.+?)\;(.+?)\;(.+?)\;(.+?)$/;
		$Domain{$domain}++;
		$Phylum{$phylum}++;
		$Class{$class}++;
		$Order{$order}++;
		$Family{$family}++;
		$Genus{$genus}++;
		$Species{$species}++;
	}
	
	my $lineage_final = "";
	my @Lineage_final = ();
	
	foreach my $domain (sort keys %Domain){
		if (($Domain{$domain} / $lineage_hit_num) >= 0.8){
			$Lineage_final[0] = $domain;
		}else{
			$Lineage_final[0] = "Not consensual";
		}
	}
		
	foreach my $phylum (sort keys %Phylum){
		if (($Phylum{$phylum} / $lineage_hit_num) >= 0.8){
			$Lineage_final[1] = $phylum;
		}else{
			$Lineage_final[1] = "Not consensual";
		}
	}

	foreach my $class (sort keys %Class){
		if (($Class{$class} / $lineage_hit_num) >= 0.8){
			$Lineage_final[2] = $class;
		}else{
			$Lineage_final[2] = "Not consensual";
		}
	}
	
	foreach my $order (sort keys %Order){
		if (($Order{$order} / $lineage_hit_num) >= 0.8){
			$Lineage_final[3] = $order;
		}else{
			$Lineage_final[3] = "Not consensual";
		}
	}	
	
	foreach my $family (sort keys %Family){
		if (($Family{$family} / $lineage_hit_num) >= 0.8){
			$Lineage_final[4] = $family;
		}else{
			$Lineage_final[4] = "Not consensual";
		}
	}	

	foreach my $genus (sort keys %Genus){
		if (($Genus{$genus} / $lineage_hit_num) >= 0.8){
			$Lineage_final[5] = $genus;
		}else{
			$Lineage_final[5] = "Not consensual";
		}
	}	

	foreach my $species (sort keys %Species){
		if (($Species{$species} / $lineage_hit_num) >= 0.8){
			$Lineage_final[6] = $species;
		}else{
			$Lineage_final[6] = "Not consensual";
		}
	}	
	
	my @Lineage_final_curated = (); # Curated final lineage
	for(my $i=0; $i<=$#Lineage_final; $i++){
		if ($Lineage_final[$i] ne "Not consensual" and $Lineage_final[$i] !~ /^\S\_\_$/){
			push @Lineage_final_curated, $Lineage_final[$i];
		}else{
			last;
		}
	}
	
	my $lineage_final_curated = join("\;",@Lineage_final_curated);
	
	return $lineage_final_curated;
}