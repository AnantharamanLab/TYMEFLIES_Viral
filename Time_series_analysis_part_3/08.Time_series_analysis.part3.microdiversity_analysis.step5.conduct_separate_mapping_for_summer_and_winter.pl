#!/usr/bin/perl

use strict;
use warnings;

# Aim: Conduct separate mapping for different seasons (summer and winter) for Fst calculations
# Note: Make sure that Bowtie2 can be called properly; Suggest to run under the environment of vRhyme (conda activate vRhyme)

# Step 1 Store IMG ID for summer and winter
my %Summer_IMG = (); # $img => 1
my %Winter_IMG = (); # $img => 1
open IN, "TYMEFLIES_metagenome_info.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $img = $tmp[0];
		my $season = $tmp[10];
		if ($season eq "Summer"){
			$Summer_IMG{$img} = 1;
		}elsif($season eq "Winter"){
			$Winter_IMG{$img} = 1;
		}
	}
}	
close IN;

# Step 2 Store all summer and winter reads address
my @Summer_reads_1 = (); # Store the address of the 1st reads
my @Summer_reads_2 = (); # Store the address of the 2nd reads
my @Winter_reads_1 = (); # Store the address of the 1st reads
my @Winter_reads_2 = (); # Store the address of the 2nd reads

open IN, "ls /storage1/Reads/TYMEFLIES_reads/*.filtered_1.fastq |";
while (<IN>){
	chomp;
	my $reads = $_;
	my ($img) = $reads =~ /reads\/(33\d+?)\.filtered/;
	if (exists $Summer_IMG{$img}){
		push @Summer_reads_1, $reads;
	}elsif (exists $Winter_IMG{$img}){
		push @Winter_reads_1, $reads;
	}
}
close IN;

open IN, "ls /storage1/Reads/TYMEFLIES_reads/*.filtered_2.fastq |";
while (<IN>){
	chomp;
	my $reads = $_;
	my ($img) = $reads =~ /reads\/(33\d+?)\.filtered/;
	if (exists $Summer_IMG{$img}){
		push @Summer_reads_2, $reads;
	}elsif (exists $Winter_IMG{$img}){
		push @Winter_reads_2, $reads;
	}
}
close IN;

# Step 3 Make mapping command
`mkdir Fst_by_inStrain_lite`;

my $ref_fasta = "All_phage_species_rep_gn_containing_AMG_n_AMG_counterpart_gene_and_flankings.fasta";
my $ref_fasta_large_idx = "All_phage_species_rep_gn_containing_AMG_n_AMG_counterpart_gene_and_flankings.bowtie2_idx";
my $threads = 5;

my $summer_reads_1 = join("\,", @Summer_reads_1);
my $summer_reads_2 = join("\,", @Summer_reads_2);
my $winter_reads_1 = join("\,", @Winter_reads_1);
my $winter_reads_2 = join("\,", @Winter_reads_2);

open OUT, ">tmp.bowtie2_mapping.sh";
print OUT "bowtie2 -x $ref_fasta_large_idx -1 $summer_reads_1 -2 $summer_reads_2 -S Fst_by_inStrain_lite/All_summer_IMG.viral_species_rep.sam -p $threads --no-unal --quiet --mm\n";
#print OUT "bowtie2 -x $ref_fasta_large_idx -1 $winter_reads_1 -2 $winter_reads_2 -S Fst_by_inStrain_lite/All_winter_IMG.viral_species_rep.sam -p $threads --no-unal --quiet --mm\n";
close OUT;

`bash tmp.bowtie2_mapping.sh`;

`rm tmp.bowtie2_mapping.sh`;

# Step 4 Convert sam file to filtered bam file
`mkdir Fst_by_inStrain_lite/original`; # Make a folder to contain the original sam and bam files

open OUT, ">tmp.covert_sam_file_to_filtered_bam_file.sh";
open IN, "ls Fst_by_inStrain_lite/All_*_IMG.viral_species_rep.sam |";
while (<IN>){
	chomp;
	my $sam = $_;
	my ($sam_name) = $sam =~ /Fst\_by\_inStrain\_lite\/(.+?)\.sam/; 
	
	# Convert sam to bam files (set the percent identity cutoff as 0.9)
	print OUT "python3 /slowdata/scripts/python_scripts/filter_coverage_file.py -s $sam -o Fst_by_inStrain_lite/$sam_name.id90.bam -p 0.90 -t 20;";
	print OUT "mv $sam Fst_by_inStrain_lite/original/$sam_name.sam; mv Fst_by_inStrain_lite/$sam_name.id90.bam Fst_by_inStrain_lite/original/$sam_name.id90.bam;";
	
	# The scaffold name of viral species representatives, which will be used as the reference for bam filtering
	my $scaf_name_file = "/storage1/data11/TYMEFLIES_phage/viral_species_representative_scaf_name.txt";
	
	# Filter bam file by using the scaffold name of only viral species representatives
	print OUT "python3 /slowdata/scripts/python_scripts/filter_bam_by_reference.py -b Fst_by_inStrain_lite/original/$sam_name.id90.bam -r $scaf_name_file -j 20;\n";
}
close OUT;

`cat tmp.covert_sam_file_to_filtered_bam_file.sh | parallel -j 2`;

`rm tmp.covert_sam_file_to_filtered_bam_file.sh`;

# Step 5 Move *viral_species_rep.id90.filtered.bam files to Fst_by_inStrain_lite folder and delete the original folder

# Step 6 Rename Fst_by_inStrain_lite to Summer_vs_Winter_Fst_analysis
