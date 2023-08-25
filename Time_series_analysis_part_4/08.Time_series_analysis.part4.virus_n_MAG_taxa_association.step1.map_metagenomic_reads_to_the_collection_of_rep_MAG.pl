#!/usr/bin/perl

use strict;
use warnings;

# Aim: Map all metagenomic reads to the collection of TYMEFLIES rep MAGs
# Note: This script should be run within conda env ViWrap-Mapping (conda activate /storage1/data11/yml_environments/ViWrap-Mapping)


# Step 1 Make batch bowtie 2 mapping command
## Step 1.1 store all IMG ID
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '^33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

## Step 1.2 Make bowtie index
#`bowtie2-build --large-index TYMEFLIES_rep_MAG.fasta TYMEFLIES_rep_MAG.bowtie2_idx --threads 20 --quiet`;

## Step 1.3 Make run command
open OUT, ">tmp.bowtie2_mapping.sh";
foreach my $img_id (sort keys %IMGID){	
	if (!( -e "$img_id/$img_id.rep_MAG.id90.bam")){
		my $reads_1 = "/storage1/Reads/TYMEFLIES_reads/$img_id.filtered_1.fastq.gz";  
		my $reads_2 = "/storage1/Reads/TYMEFLIES_reads/$img_id.filtered_2.fastq.gz";  
		my $threads = 1;
	
		print OUT "bowtie2 -x TYMEFLIES_rep_MAG.bowtie2_idx -1 $reads_1 -2 $reads_2 -S $img_id/$img_id.rep_MAG.sam -p $threads --no-unal --quiet --mm;";
		print OUT "samtools view -bS $img_id/$img_id.rep_MAG.sam > $img_id/$img_id.rep_MAG.bam -@ $threads 2> /dev/null;";
		print OUT "samtools flagstat $img_id/$img_id.rep_MAG.bam > $img_id/$img_id.rep_MAG.stat 2> /dev/null;";
		print OUT "python3 /slowdata/scripts/python_scripts/filter_coverage_file.py -b $img_id/$img_id.rep_MAG.bam -o $img_id/$img_id.rep_MAG.id90.bam -p 0.90 -t $threads;";
		print OUT "samtools flagstat $img_id/$img_id.rep_MAG.id90.bam > $img_id/$img_id.rep_MAG.id90.stat 2> /dev/null;";
		print OUT "samtools index $img_id/$img_id.rep_MAG.id90.bam 2> /dev/null;";
		print OUT "rm $img_id/$img_id.rep_MAG.sam $img_id/$img_id.rep_MAG.sorted.bam $img_id/$img_id.rep_MAG.bam\n";
	}
}
close OUT;

#`cat tmp.bowtie2_mapping.sh | parallel -j 30`;

#`rm tmp.bowtie2_mapping.sh`;

