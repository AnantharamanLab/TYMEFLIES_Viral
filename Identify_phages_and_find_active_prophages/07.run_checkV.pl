#!/usr/bin/perl

use strict;
use warnings;

# this script should be run within conda env "CheckV_v0.8.1"
# before using, do this command to indicate database: export CHECKVDB=/slowdata/databases/checkv-db-v0.6

my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

open OUT, ">tmp.run_checkV.sh";
foreach my $img_id (sort keys %IMGID){
	if (!(-d "$img_id/CheckV_phage_scaffold")){
		print OUT "checkv end_to_end $img_id/VIBRANT_$img_id.a.v2.min2000/$img_id.a.phage_scaffold.fna $img_id/CheckV_phage_scaffold -t 1\n";
	}
}
close OUT;

`cat tmp.run_checkV.sh | parallel -j 15`;

`rm tmp.run_checkV.sh`;
