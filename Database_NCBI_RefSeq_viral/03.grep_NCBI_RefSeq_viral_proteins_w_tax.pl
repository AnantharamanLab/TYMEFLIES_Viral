#!/usr/bin/perl

use strict;
use warnings;

# AIM: Grep NCBI RefSeq viral proteins with tax information parsed

# Step 1. Store the NCBI_RefSeq_viral_pro_w_tax_list
my %NCBI_RefSeq_viral_pro_w_tax_list = ();
open IN, "viral.protein.tax.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	$NCBI_RefSeq_viral_pro_w_tax_list{$tmp[0]} = 1;
}
close IN;

# Step 2. Store all the protein sequences and change the head 
my %Seq = _store_seq("viral.protein.faa");

my %Seq2 = (); # Change the head: >YP_009137152.1 => YP_009137152
foreach my $key (sort keys %Seq){
	my ($key_new) = $key =~ /^(>.+?)\./;
	$Seq2{$key_new} = $Seq{$key};
}

# Step 3. Only keep proteins with tax information parsed
foreach my $key (sort keys %Seq2){
	my ($key_clean) = $key =~ /^>(.+?)$/;
	if (! exists $NCBI_RefSeq_viral_pro_w_tax_list{$key_clean}){
		delete $Seq2{$key};
	}
}

# Step 4. Write down the result
open OUT, ">viral.protein.w_tax.faa";
foreach my $key (sort keys %Seq2){
	print OUT "$key\n$Seq2{$key}\n";
}
close OUT;
# subroutine

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