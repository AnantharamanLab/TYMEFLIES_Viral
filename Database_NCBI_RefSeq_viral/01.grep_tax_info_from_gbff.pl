#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

# AIM: To get taxonomy information from input gb file

my $in  = Bio::SeqIO->new(-file => "viral.genomic.gbff",
                         -format => 'genbank');

open OUT, ">viral.genomic.tax.txt";
while (my $seq = $in->next_seq()){
	my @Full_tax = $seq->species->classification;
	@Full_tax = reverse(@Full_tax);
	my $species = $seq->species->node_name;
	my $full_tax = join("\;",@Full_tax);
	my $id = $seq->id;
	print OUT "$id\t$full_tax\t$species\n";	
}
close OUT;

