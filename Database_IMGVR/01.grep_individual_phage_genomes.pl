#!/usr/bin/perl

use strict;
use warnings;

#########################################################
# Step 1. Check phage genome name congruency            #
# Step 2. Make individual phage genomes in a new folder #
#########################################################

## Step 1. Check phage genome name congruency      

# Step 1.1 Determine the genome name from IMGVR_all_nucleotides_header file
my %Gn2scf = (); # $gn => $scf separated by "\t"
open IN, "IMGVR_all_nucleotides_header.txt";
while (<IN>){
	chomp;
	my $scf = $_;
	my $gn = "";
	
	# Situation 1. Scf contains UViG
	# All scf has 2 or 3 "|", the content before the 1st "|" is the genome name
	if ($scf =~ /UViG/){
		($gn) = $scf =~ /^(.+?)\|/;
	}else{
		# Sitation 2. Scf has IMG ID inside, and before IMG ID, the content is not IMG ID again, the content is the genome name
		if ($scf =~ /^.+?\|\d\d\d\d\d\d\d\d\d\d\|/ and $scf !~ /^\d\d\d\d\d\d\d\d\d\d\|\d\d\d\d\d\d\d\d\d\d\|/){
			($gn) = $scf =~ /^(.+?)\|\d\d\d\d\d\d\d\d\d\d\|/;
		}elsif (/^\d\d\d\d\d\d\d\d\d\d\|/){ # Sitation 3. Scf has IMG ID inside and in the front, and it is a single-scaffold phage name, the phage genome name is the same as the scaffold name
			$gn = $scf;
		}elsif(!/\|/){ # Situation 4. Scf doesn't have "|" inside, the whole scf name is the phage genome name
			$gn = $scf;
		}else{ # Other situations, manually check
			$gn = $scf;
		}
	}

	if ($gn){
		if (!exists $Gn2scf{$gn}){
			$Gn2scf{$gn} = $scf;
		}else{
			$Gn2scf{$gn} .= "\t".$scf;
		}
	}

}
close IN;

# Step 1.2 Store the genome name from "IMGVR_all_Sequence_information.tsv", and check if our genome name is the same as this
my %Gn_from_all_seq_info = (); # $gn => 1
open IN, "IMGVR_all_Sequence_information.tsv";
while (<IN>){
	chomp;
	if (!/^#/){
		my @tmp = split (/\t/);
		my $gn = $tmp[0];
		$Gn_from_all_seq_info{$gn} = 1;
	}
}
close IN;

# Check if our gn ids are all present in theirs
foreach my $gn1 (sort keys %Gn2scf){
	if (!exists $Gn_from_all_seq_info{$gn1}){
		print "$gn1 is not in IMGVR_all_Sequence_information.tsv\n";
	}
}

# Check if their gn ids are all present in ours
foreach my $gn2 (sort keys %Gn_from_all_seq_info){
	if (!exists $Gn2scf{$gn2}){
		print "$gn2 is not in our Gn2scf hash\n";
	}
}

# Check the numbers of elements in both hash
my $element_num_Gn2scf = scalar (keys %Gn2scf);
my $element_num_Gn_from_all_seq_info = scalar (keys %Gn_from_all_seq_info);
print "Our hash contains $element_num_Gn2scf elements\n";
print "Their hash contains $element_num_Gn_from_all_seq_info elements\n";

## Step 2. Make individual phage genomes in a new folder

`mkdir Individual_phage_genomes`;
my %Seq = _store_seq("IMGVR_all_nucleotides.fna");

foreach my $gn (sort keys %Gn2scf){
	my %Each_fasta_seq = (); # The seq hash for each fasta file
	
	my @Scf = split (/\t/, $Gn2scf{$gn});
	foreach my $scf (@Scf){
		$Each_fasta_seq{">$scf"} = $Seq{">$scf"};
	}
	
	# Replace "|" to "__" in the gn name
	my $gn_new = $gn; $gn_new =~ s/\|/\_\_/g;
	open OUT, ">Individual_phage_genomes/$gn_new.fasta";
	foreach my $key (sort keys %Each_fasta_seq){
		# Replace "|" to "__" in sequence header
		my $key_new = $key; $key_new =~ s/\|/\_\_/g;
		print OUT "$key_new\n$Each_fasta_seq{$key}\n";
	}
	close OUT;
}

# Step 2.1 write down "IMGVR_all_genome_scaffolds.txt"
open OUT, ">IMGVR_all_genome_scaffolds.txt";
print OUT "viral genome\tscaffold number\tscaffolds\n";
foreach my $gn (sort keys %Gn2scf){
	my $scf_num = 0;
	my @Scf = split (/\t/, $Gn2scf{$gn});
	$scf_num = scalar @Scf;
	print OUT "$gn\t$scf_num\t$Gn2scf{$gn}\n";
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

