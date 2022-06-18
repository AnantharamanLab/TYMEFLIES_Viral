#!/usr/bin/perl

use strict;
use warnings;

# Aim: Map all metagenomic reads to the collection of species representatives
# Note: This script should be run within conda env vRhyme (conda activate vRhyme)

# Step 1 Get the viral genomic sequences of all species representatives
my %All_viral_genome_seq = _store_seq("/storage1/data11/TYMEFLIES_phage/Host_prediction/All_phage_genomes.fasta");

## Step 1.1 Store species info
my %Species = (); # $gn_rep => $gns 
open IN, "/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt";
while (<IN>){
	chomp;
	my @tmp = split (/\t/);
	my $gn_rep = $tmp[0];
	my $gns = $tmp[1];
	$Species{$gn_rep} = $gns;
}
close IN;

## Step 1.2 Store AMG KO information
my %AMG_summary = (); # $pro => $ko
my %KOs= (); # $ko => 1;
my %IMG2date = (); # $img_id => $date_n_season
my %AMG_containing_viral_gn = (); # $gn => 1
open IN, "AMG_analysis/AMG_summary.txt";
while (<IN>){
	chomp;
	if (!/^Pro/){
		my @tmp = split (/\t/);
		my $pro = $tmp[0];
		my $ko = $tmp[2];
		my $ko_detail = $tmp[3];
		my $date_n_season = $tmp[1];
		$AMG_summary{$pro} = $ko;
		my ($img_id) = $pro =~ /^(33.+?)\_/;
		$IMG2date{$img_id} = $date_n_season;
		$KOs{$ko} = $ko_detail;
		
		my ($gn) = $pro =~ /^(.+?\_\_.+?)\_\_/;
		$AMG_containing_viral_gn{$gn} = 1;
	}
}
close IN;
=pod
open OUT, ">All_phage_species_rep_gn_containing_AMG.fasta";
foreach my $header (sort keys %All_viral_genome_seq){
	my ($gn) = $header =~ /^>(.+?\_\_.+?)\_\_/;
	if (exists $Species{$gn} and exists $AMG_containing_viral_gn{$gn}){
		print OUT "$header\n$All_viral_genome_seq{$header}\n";
	}
}
close OUT;
=cut
# Step 2 Make batch bowtie 2 mapping command
## Step 2.1 store all IMG ID
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '^33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

## Step 2.2 Make bowtie index
#`cat AMG_counterpart_genes_and_flankings.fasta reference_fasta_for_metapop/All_phage_species_rep_gn_containing_AMG.fasta > All_phage_species_rep_gn_containing_AMG_n_AMG_counterpart_gene_and_flankings.fasta`;
#`bowtie2-build --large-index All_phage_species_rep_gn_containing_AMG_n_AMG_counterpart_gene_and_flankings.fasta All_phage_species_rep_gn_containing_AMG_n_AMG_counterpart_gene_and_flankings.bowtie2_idx --threads 10 --quiet`;

## Step 2.3 Make run command
open OUT, ">tmp.bowtie2_mapping.sh";
foreach my $img_id (sort keys %IMGID){	
	if (!( -e "$img_id/$img_id.viral_species_rep.id90.bam")){
		my $reads_1 = "/storage1/Reads/TYMEFLIES_reads/$img_id.filtered_1.fastq";  
		my $reads_2 = "/storage1/Reads/TYMEFLIES_reads/$img_id.filtered_2.fastq";  
		my $threads = 1;
	
		print OUT "bowtie2 -x All_phage_species_rep_gn_containing_AMG_n_AMG_counterpart_gene_and_flankings.bowtie2_idx -1 $reads_1 -2 $reads_2 -S $img_id/$img_id.viral_species_rep.sam -p $threads --no-unal --quiet --mm;";
		print OUT "samtools view -bS $img_id/$img_id.viral_species_rep.sam > $img_id/$img_id.viral_species_rep.bam -@ $threads 2> /dev/null;";
		print OUT "samtools flagstat $img_id/$img_id.viral_species_rep.bam > $img_id/$img_id.viral_species_rep.stat 2> /dev/null;";
		print OUT "python3 /slowdata/scripts/python_scripts/filter_coverage_file.py -b $img_id/$img_id.viral_species_rep.bam -o $img_id/$img_id.viral_species_rep.id90.bam -p 0.90 -t $threads;";
		print OUT "samtools flagstat $img_id/$img_id.viral_species_rep.id90.bam > $img_id/$img_id.viral_species_rep.id90.stat 2> /dev/null;";
		#print OUT "coverm contig --methods metabat --bam-files $img_id/$img_id.viral_species_rep.id90.bam > $img_id/$img_id.viral_species_rep.id90.coverm_depth.txt 2> /dev/null;";
		print OUT "samtools index $img_id/$img_id.viral_species_rep.id90.bam 2> /dev/null;";
		print OUT "rm $img_id/$img_id.viral_species_rep.sam $img_id/$img_id.viral_species_rep.sorted.bam $img_id/$img_id.viral_species_rep.bam\n";
	}
}
close OUT;

`cat tmp.bowtie2_mapping.sh | parallel -j 40`;

`rm tmp.bowtie2_mapping.sh`;



## Subroutine

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