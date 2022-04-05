#!/usr/bin/perl

use strict;
use warnings;

# Aim: mapping reads of each sample onto all scaffolds within this sample
# Note: This should be run under the vRhyme conda env by "conda activate vRhyme"

# 1. store all IMG ID
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '^33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

# 2. bowtie2 mapping 
# 2.1 store the metagenome info
my %IMG2info = (); 
open IN, "TYMEFLIES_metagenome_info.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my $line = $_;
		my @tmp = split (/\t/);
		my $img_id = $tmp[0];
		$IMG2info{$img_id} = $line; 
	}
}
close IN;

# 2.2 run bowtie2 for all samples
open OUT, ">tmp.bowtie2_mapping.sh";
foreach my $img_id (sort keys %IMGID){	
	if (!( -e "$img_id/$img_id.id97.bam")){
		my $reads_1 = "/storage1/Reads/TYMEFLIES_reads/$img_id.filtered_1.fastq";  
		my $reads_2 = "/storage1/Reads/TYMEFLIES_reads/$img_id.filtered_2.fastq";  
		my $threads = 1;
		my $fna = "$img_id/$img_id.a.fna";
		
		print OUT "bowtie2-build $fna $img_id/$img_id.bowtie2_idx --quiet;";	
		print OUT "bowtie2 -x $img_id/$img_id.bowtie2_idx -1 $reads_1 -2 $reads_2 -S $img_id/$img_id.sam -p $threads --no-unal --quiet;";
		print OUT "python3 /slowdata/scripts/python_scripts/filter_coverage_file.py -s $img_id/$img_id.sam -o $img_id/$img_id.id97.bam -t 1;";
		print OUT "coverm contig --methods metabat --bam-files $img_id/$img_id.id97.bam > $img_id/$img_id.id97.coverm_depth.txt 2> /dev/null;";
		print OUT "samtools index $img_id/$img_id.id97.bam 2> /dev/null;";
		print OUT "rm $img_id/$img_id.sam $img_id.bam $img_id.sorted.bam\n";
	}
}
close OUT;

my $threads = 50;

`cat tmp.bowtie2_mapping.sh | parallel -j $threads`;

`rm tmp.bowtie2_mapping.sh`;





