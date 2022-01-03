#!/usr/bin/perl

use strict;
use warnings;

# this script should be run under the conda env "BBTools_v37.62"

my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

# print bbmap mapping scripts
#open OUT, ">tmp.bbmap_mapping.sh";
foreach my $img_id (sort keys %IMGID){
	if (! (-e "$img_id\/$img_id\.a\.depth\.txt")){
	
		my $fna = "$img_id/$img_id.a.fna";
		my $reads = "/storage1/Reads/TYMEFLIES_reads/$img_id.filtered.fq.gz";
		my $threads = 20;
		
		`bbmap.sh ref=$fna in=$reads outm=$img_id.m.sam nodisk=f interleaved=t ambiguous=random covstats=$img_id.covstat twocolumn=t t=$threads; rm -r ref; rm $img_id.m.sam`;
		
	}
}
#close OUT;





