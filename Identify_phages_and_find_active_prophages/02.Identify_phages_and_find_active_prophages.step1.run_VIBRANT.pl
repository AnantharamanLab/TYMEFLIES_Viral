#!/usr/bin/perl

use strict;
use warnings;

# this script should be run under 'vibrant_env' conda env

my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

# print VIBRANT scripts
open OUT, ">tmp.run_VIBRANT.sh";
foreach my $img_id (sort keys %IMGID){
	if (! (-d "$img_id/VIBRANT\_${img_id}\.a")){
		my $fna = "$img_id/$img_id.a.fna";
		print OUT "python3 /home/zhichao/miniconda3/envs/vibrant_env/bin/VIBRANT_run.py -i $img_id/$img_id.a.fna -folder $img_id -t 40 -d /slowdata/data4/VIBRANT/VIBRANT_v1.2.1/databases -m /slowdata/data4/VIBRANT/VIBRANT_v1.2.1/files\n";
	}
}
close OUT;

`bash tmp.run_VIBRANT.sh`;

`rm tmp.run_VIBRANT.sh`;
