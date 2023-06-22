#!/usr/bin/perl

use strict;
use warnings;

# this script should be run within conda env "vrhyme" by run "conda activate /home/kieft/miniconda3/envs/vrhyme" or "conda activate vrhyme"

my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

# run vRhyme
open OUTT, ">tmp.run_vRhyme.sh";
foreach my $img_id (sort keys %IMGID){
	#$img_id = "3300020571";
	# test if the result folder is present
	if (!(-d "$img_id/vRhyme_result") and -d "$img_id\/VIBRANT\_$img_id\.a\/VIBRANT\_results\_$img_id.a"){
		# deinterleave reads
		my $read_dir = "/storage1/Reads/TYMEFLIES_reads";
				
		# run vRhyme
		my $threads = 1;
		print OUTT "vRhyme -i $img_id/VIBRANT_$img_id.a.v2.min2000/$img_id.a.phage_scaffold.fna -t $threads -r $read_dir/$img_id.filtered*.fastq.gz -o $img_id/vRhyme_result\n";
	}
}
close OUTT;

`cat tmp.run_vRhyme.sh | parallel -j 20`;

`rm tmp.run_vRhyme.sh`;