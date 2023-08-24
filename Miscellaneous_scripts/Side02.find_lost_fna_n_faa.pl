#!/usr/bin/perl

use strict;
use warnings;

my %MAG2IMGID = (); #$mag => $img_id
my %IMGID = (); #$img_id => $mag separated by "\t"
open IN, "ls *.tar.gz |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($mag) = $file =~ /^(.+?)\.tar\.gz/;
	my ($img_id) = $mag =~ /^(.+?)\_/;
	$MAG2IMGID{$mag} = $img_id;
	if (exists $IMGID{$img_id}){
		$IMGID{$img_id} .= "\t".$mag;
	}else{
		$IMGID{$img_id} = $mag;
	}
}
close IN;

# print lost files
my %Lost = (); # $mag => $file separated by "\t"
open IN, "ls  -l | grep ^d| grep '33' | sed -r s'/ +/\t/g' | cut -f 9 |";
while (<IN>){
	chomp;
	my $img_id = $_;
	my @MAGs = split (/\t/, $IMGID{$img_id});
	foreach my $mag (@MAGs){
		my $fna = "$img_id\/$mag\.fna";
		my $faa = "$img_id\/$mag\.faa";
#		my $gff = "$img_id\/$mag\.gff";
		my $cog = "$img_id\/$mag\.cog\.txt";
		my $pfam = "$img_id\/$mag\.pfam\.txt";
		my $phylodist = "$img_id\/$mag\.phylodist\.txt";
		my $ko = "$img_id\/$mag\.ko\.txt";
		my $ec = "$img_id\/$mag\.ec\.txt";
		my $gene_product = "$img_id\/$mag\.gene\_product\.txt";
		my $mbin = "$img_id\/mbin\_datafile\_${img_id}\.txt";
		
		if (! (-e "$fna")){
			if (!exists $Lost{$mag}){
				$Lost{$mag} = "fna";
			}else{
				$Lost{$mag} .= "\tfna";
			}
		}
		
		if (! (-e "$faa")){
			if (!exists $Lost{$mag}){
				$Lost{$mag} = "faa";
			}else{
				$Lost{$mag} .= "\tfaa";
			}
		}
=pod		
		if (! (-e "$gff")){
			if (!exists $Lost{$mag}){
				$Lost{$mag} = "gff";
			}else{
				$Lost{$mag} .= "\tgff";
			}
		}
=cut
		if (! (-e "$cog")){
			if (!exists $Lost{$mag}){
				$Lost{$mag} = "cog";
			}else{
				$Lost{$mag} .= "\tcog";
			}
		}	

		if (! (-e "$pfam")){
			if (!exists $Lost{$mag}){
				$Lost{$mag} = "pfam";
			}else{
				$Lost{$mag} .= "\tpfam";
			}
		}	

		if (! (-e "$phylodist")){
			if (!exists $Lost{$mag}){
				$Lost{$mag} = "phylodist";
			}else{
				$Lost{$mag} .= "\tphylodist";
			}
		}	

		if (! (-e "$ko")){
			if (!exists $Lost{$mag}){
				$Lost{$mag} = "ko";
			}else{
				$Lost{$mag} .= "\tko";
			}
		}	

		if (! (-e "$ec")){
			if (!exists $Lost{$mag}){
				$Lost{$mag} = "ec";
			}else{
				$Lost{$mag} .= "\tec";
			}
		}	

		if (! (-e "$gene_product")){
			if (!exists $Lost{$mag}){
				$Lost{$mag} = "gene\_product";
			}else{
				$Lost{$mag} .= "\tgene\_product";
			}
		}	

		if (! (-e "$mbin")){
			if (!exists $Lost{$mag}){
				$Lost{$mag} = "mbin";
			}else{
				$Lost{$mag} .= "\tmbin";
			}
		}		
	}
}
close IN;


foreach my $mag (sort keys %Lost){
	print "$mag\t$Lost{$mag}\n";
	#`cp $mag.tar.gz tmp`;
}








=pod
<bin_oid>.fna - FASTA nucleic acid file of scaffolds for genome bin.
<bin_oid>.faa - FASTA amino acid file for genome bin.
<bin_oid>.cog.txt - Tab delimited file for COG annotation.
<bin_oid>.pfam.txt - Tab delimited file for Pfam annotation.
<bin_oid>.phylodist.txt - Tab delimited file for phylo distribution
<bin_oid>.ko.txt - Tab delimited file for KO annotation.
<bin_oid>.ec.txt - Tab delimited file for EC annotation.
<bin_oid>.gene_product.txt - Product name assignment.
=cut