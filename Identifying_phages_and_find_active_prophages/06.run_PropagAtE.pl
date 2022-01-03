#!/usr/bin/perl

use strict;
use warnings;
use File::Slurp;

# this script should be run within conda env "vRhyme"; vRhyme env includes all the dependencies of PropogAtE

my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

# run 1 metagenome to have a test
#%IMGID = (); $IMGID{"3300034109"} = 1;

# run propagate
open OUT, ">tmp.run_propagate.sh";
foreach my $img_id (sort keys %IMGID){
	# test if prophage is present
	if (-d "$img_id\/VIBRANT\_$img_id\.a" and -e "$img_id/VIBRANT_$img_id.a.v2.min2000/$img_id.a.prophage_scaffold.fna" and !( -d "$img_id\/PropagAtE\_result")){
			
			# deinterleave reads
			my $read_dir = "/storage1/Reads/TYMEFLIES_reads";
			#`gzip -dc $read_dir/$img_id.filtered.fastq.gz | bash deinterleave_fastq.sh $read_dir/$img_id.filtered.1.fastq $read_dir/$img_id.filtered.2.fastq`;
			
			# prepare prophage scaffold fasta file
            my %Fna = _store_seq("$img_id/$img_id.a.fna");
            my %Prophage_scf_seq = ();
            open IN, "$img_id/VIBRANT_$img_id.a.v2.min2000/prophage_coordinates.txt";
            while (<IN>){
                chomp;
                if (!/^scaffold/){
                   my @tmp = split (/\t/);
                   my $scf = $tmp[0]; $scf = ">".$scf;
                   $Prophage_scf_seq{$scf} = $Fna{$scf};
                }
            }
            close IN;

			if (%Prophage_scf_seq){
				open OUT2, ">$img_id/Prophage_scf.fna";
				foreach my $scf (sort keys %Prophage_scf_seq){
					print OUT2 "$scf\n$Prophage_scf_seq{$scf}\n";
				}
				close OUT2;
			}
			
			# run propagate
			my $threads = 1;
			if (-e "$img_id/Prophage_scf.fna"){
				print OUT "/slowdata/data4/PropagAtE_prophage_activity/PropagAtE_GitHub_v1.1.0/Propagate/Propagate -f $img_id/Prophage_scf.fna -r $read_dir/$img_id.filtered_1.fastq $read_dir/$img_id.filtered_2.fastq -v $img_id/VIBRANT_$img_id.a.v2.min2000/prophage_coordinates.txt -o $img_id/PropagAtE_result --clean -t $threads\n";
			}
	
	}
}
close OUT;

`cat tmp.run_propagate.sh | parallel -j 20`;

`rm tmp.run_propagate.sh`;

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





