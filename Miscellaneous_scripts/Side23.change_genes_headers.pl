#!/usr/bin/perl
  
use strict;
use warnings;

# Aim: Change the gene file headers to run MetaPop
# Three places to change: (1) The start and stop positions of each gene from prophages were changed
#                         (2) The stop position of the last gene of prophage should be (the length of prophage - 1), 
#							  and the last nucleotide of the sequence of last gene is also deleted
#						  (3) The protein ID start from 1 for prophage genes
#						  (4) The partial and start_type in the headers are changed
#						  (5) The gene IDs in the header for all genes are changed
#						  (6) The order of genes are sorted according to the gene ID
 
my %Seq = _store_seq("/storage1/data11/TYMEFLIES_phage/All_phage_species_rep_gn.genes");

# Step 1 Store the VIBRANT prophage coordinates
my %Prophage_fragment2start_stop = (); # $prophage_fragment (i.e., "Ga0453920_0001995_fragment_1") => $start."\_".$stop."\_".$protein_start
open IN, "/storage1/data11/TYMEFLIES_phage/VIBRANT_prophage_coordinates.txt";
while (<IN>){
	chomp;
	if (!/^scaffold/){
		my @tmp = split (/\t/);
		my $prophage_fragment = $tmp[1];
		my $protein_start = $tmp[2];
		my $start = $tmp[5];
		my $stop = $tmp[6];
		$Prophage_fragment2start_stop{$prophage_fragment} = $start."\_".$stop."\_".$protein_start;
	}
}
close IN;

# Step 2 Parse gene file
my %Seq_new = %Seq; # Store a new Seq hash
my %New_gene2old_gene_map = (); # $gene_new => $gene
foreach my $header (sort keys %Seq){
	my @tmp = split (/\s\#\s/, $header);
	
	my ($gene) = $tmp[0] =~ /^>(.+?)$/;
	my $start = $tmp[1]; my $stop = $tmp[2];
	my $direction = $tmp[3]; my $rest = $tmp[4];
	
	if ($gene =~ /fragment/){
		my ($prophage_fragment) = $gene =~ /(Ga.+?fragment\_.+?)\_/;
		my ($start_at_whole_scaffold, $stop_at_whole_scaffold, $protein_start) = $Prophage_fragment2start_stop{$prophage_fragment} =~ /^(.+?)\_(.+?)\_(.+?)$/;
		my $length_at_whole_scaffold = $stop_at_whole_scaffold - $start_at_whole_scaffold + 1;
		
		# Change the start and stop position of each gene from prophages
		$start = $start - ($start_at_whole_scaffold - 1);
		$stop = $stop - ($start_at_whole_scaffold - 1);
		
		# The stop position of the last gene of prophage should be (the length of prophage - 1)
		if ($stop == $length_at_whole_scaffold){
			$stop = $length_at_whole_scaffold - 1;
		}
		
		# Change the start protein ID to 1 for prophage genes
		my ($scaffold, $protein_id) = $gene =~ /^(.+)\_(\d+?)$/;
		$protein_id = $protein_id - ($protein_start - 1);
		my $gene_new = $scaffold."\_".$protein_id;
		$New_gene2old_gene_map{$gene_new} = $gene;
		
		# Make the new header
		my $header_new = ">".$gene_new." # ".$start." # ".$stop." # ".$direction." # ".$rest;
		
		# Store the header and sequence into the new hash
		delete $Seq_new{$header};
		
		if ($stop == $length_at_whole_scaffold - 1){
			if ($direction eq "1"){
				my $seq = $Seq{$header};
				# Delete the last nucleotide
				my $seq_length = length($seq);
				$seq = substr($seq, 0, ($seq_length - 1));
				
				# Change the partial and start_type of header_new
				$header_new =~ s/partial\=00/partial\=01/g; # Change partial
								
				$Seq_new{$header_new} = $seq;
			}elsif($direction eq "-1"){
				my $seq = $Seq{$header};
				# Delete the first nucleotide
				my $seq_length = length($seq);
				$seq = substr($seq, 1, ($seq_length - 1));
				
				# Change the partial and start_type of header_new
				$header_new =~ s/partial\=00/partial\=01/g; # Change partial
				$header_new =~ s/start\_type\=.TG/start\_type\=Edge/g; # Change start_type
				
				$Seq_new{$header_new} = $seq;
			}
		}else{
			$Seq_new{$header_new} = $Seq{$header};
		}
	}
}

## Change the gene ID for all genes (including both prophage and other non-prophage genes)
my %Seq_new_2 = ();
foreach my $header_new (sort keys %Seq_new){
	my $seq_new = $Seq_new{$header_new};
	my ($gene_id) = $header_new =~ /\_\_(Ga.+?)\s\#/;
	$header_new =~ s/ID\=.+?\;/ID\=$gene_id\;/g;
	$Seq_new_2{$header_new} = $seq_new;
}

# Step 3 Write down the sequences
my @Sorted_Seq_new_2 = sort { # Sort the header
    my ($a1, $a2) = $a =~ /^>(.+?)\_(\d+?)\s\#\s/;
    my ($b1, $b2) = $b =~ /^>(.+?)\_(\d+?)\s\#\s/;
	$a2 = (sprintf "%04d", $a2);
	$b2 = (sprintf "%04d", $b2);
	my $aa = $a1."\_".$a2; my $bb = $b1."\_".$b2; 
    $aa cmp $bb;
} keys %Seq_new_2;

open OUT, ">/storage1/data11/TYMEFLIES_phage/All_phage_species_rep_gn.mdfed.genes";
foreach my $key (@Sorted_Seq_new_2){
	print OUT "$key\n$Seq_new_2{$key}\n";
}
close OUT;

# Step 4 Write down the new gene to old gene map
open OUT, ">/storage1/data11/TYMEFLIES_phage/New_gene2old_gene_map.txt";
foreach my $key (sort keys %New_gene2old_gene_map){
	print OUT "$key\t$New_gene2old_gene_map{$key}\n";
}
close OUT;



## Subroutine

sub _store_seq{  # Store the full header
	my $file = $_[0];
	my %Seq = (); my $head = "";
	open _IN, "$file";
	while (<_IN>){
		chomp;
		if (/>/){
			($head) = $_ =~ /^(>.+?)$/;
			$Seq{$head} = "";	
		}else{
			$Seq{$head} .= $_;
		}
	}
	close _IN;
	return %Seq;
}
