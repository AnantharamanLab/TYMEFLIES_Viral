#!/usr/bin/perl

use strict;
use warnings;

# AIM: Run blastn to filter viral sequences from MAGs (both TYMEFLIES and GEM MAGs)
=pod
# Step 1. concatenate all nt sequences from both databases and split them into 20000 seq-containing fsa files
my $work_dir_TYMEFLIES_MAGs = "/storage1/data11/TYMEFLIES_phage/Binning_Data";
`find $work_dir_TYMEFLIES_MAGs/33* -name '*.fna' -exec cat {} + > $work_dir_TYMEFLIES_MAGs/TYMEFLIES_MAGs_all_seq.fna`;
my $work_dir_GEM_MAGs = "/slowdata/databases/GEM";
`find  $work_dir_GEM_MAGs/fna -name '*.fna' -exec cat {} + > $work_dir_GEM_MAGs/GEM_MAGs_all_seq.fna`;

`mkdir $work_dir_TYMEFLIES_MAGs/split_fsa`;
`perl /storage1/data11/TYMEFLIES_phage/split_multifasta.pl --in $work_dir_TYMEFLIES_MAGs/TYMEFLIES_MAGs_all_seq.fna --output_dir=$work_dir_TYMEFLIES_MAGs/split_fsa --seqs_per_file=20000`;

`mkdir $work_dir_GEM_MAGs/split_fsa`;
`perl /storage1/data11/TYMEFLIES_phage/split_multifasta.pl --in $work_dir_GEM_MAGs/GEM_MAGs_all_seq.fna --output_dir=$work_dir_GEM_MAGs/split_fsa --seqs_per_file=20000`;

`rm $work_dir_TYMEFLIES_MAGs/TYMEFLIES_MAGs_all_seq.fna $work_dir_GEM_MAGs/GEM_MAGs_all_seq.fna`;


# Step 2. Make viral sequence databases (blastn database)
# Viral sequence databases include: 1) IMG VR v3 all phages 2) NCBI RefSeq all viruses
#                                   3) TYMEFLIES all phages 4) Lake Mendota time series (2008-2012) all phages 

`mkdir All_viral_seq_db`;

#/slowdata/databases/IMGVR-NCBI_phages/IMGVR_all_nucleotides.fna  -> 2377994 sequences
#/slowdata/databases/NCBI_RefSeq_viral/viral.genomic.fna -> 14724 sequences
`find /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed/ -name '*.fasta' -exec cat {} + > All_viral_seq_db/TYMEFLIES_all_phages.fasta`; # 1804721 sequences
`find /storage1/data11/LakeMendota_2008_to_2012/33*/vRhyme_best_bins_fasta_parsed/ -name '*.fasta' -exec cat {} + > All_viral_seq_db/Lake_Mendota_2008_to_2012_all_phages.fasta`; # 34071 sequences

`cat /slowdata/databases/IMGVR-NCBI_phages/IMGVR_all_nucleotides.fna /slowdata/databases/NCBI_RefSeq_viral/viral.genomic.fna All_viral_seq_db/TYMEFLIES_all_phages.fasta All_viral_seq_db/Lake_Mendota_2008_to_2012_all_phages.fasta > All_viral_seq_db/All_viral_seq.fasta`;
# 4231510 sequences in All_viral_seq_db/All_viral_seq.fasta

`rm All_viral_seq_db/TYMEFLIES_all_phages.fasta All_viral_seq_db/Lake_Mendota_2008_to_2012_all_phages.fasta`;

`makeblastdb -in All_viral_seq_db/All_viral_seq.fasta -title All_viral_seq -dbtype nucl -out All_viral_seq_db/All_viral_seq_blastdb`;
=cut
# Step 3. Run blastn 
# Step 3.1 Write down blastn command for TYMEFLIES MAGs and GEM MAGs
`mkdir Host_prediction`;
`mkdir Host_prediction/filter_MAGs_out`;

open OUT, ">tmp.filter_MAGs_blastn_1.sh";
open IN, "ls /storage1/data11/TYMEFLIES_phage/Binning_Data/split_fsa/*.fsa |";
while (<IN>){
	chomp;
	my $fsa_file = $_;
	my ($out_name) = $fsa_file =~ /split_fsa\/(.+?)\.fsa/;
	$out_name = "TYMEFLIES_MAGs_".$out_name;
	if (!(-e "Host_prediction/filter_MAGs_out/$out_name.blastn_out.txt")){
		print OUT "blastn -task megablast -evalue 0.001 -max_target_seqs 25000 -perc_identity 90 -num_threads 1 -outfmt 6 -query $fsa_file -db All_viral_seq_db/All_viral_seq_blastdb -out Host_prediction/filter_MAGs_out/$out_name.blastn_out.txt\n";
	}
}
close IN;
close OUT;


open OUT, ">tmp.filter_MAGs_blastn_2.sh";
open IN, "ls /slowdata/databases/GEM/split_fsa/*.fsa |";
while (<IN>){
	chomp;
	my $fsa_file = $_;
	my ($out_name) = $fsa_file =~ /split_fsa\/(.+?)\.fsa/;
	$out_name = "GEM_MAGs_".$out_name;
	if (!(-e "Host_prediction/filter_MAGs_out/$out_name.blastn_out.txt")){
		print OUT "blastn -task megablast -evalue 0.001 -max_target_seqs 25000 -perc_identity 90 -num_threads 1 -outfmt 6 -query $fsa_file -db All_viral_seq_db/All_viral_seq_blastdb -out Host_prediction/filter_MAGs_out/$out_name.blastn_out.txt\n";
	}
}
close IN;
close OUT;

# Step 3.2 Run the command with 5 cpus

`cat tmp.filter_MAGs_blastn_1.sh | parallel -j 10`;
`cat tmp.filter_MAGs_blastn_2.sh | parallel -j 10`;

`rm tmp.filter_MAGs_blastn_1.sh`;
`rm tmp.filter_MAGs_blastn_2.sh`;
=pod
# Step 4. Parse result and provide filtered Genome and filtered contigs list
## Step 4.1 Store all contig to MAG info
my %TYMEFLIES_contig2MAG = ();
my %TYMEFLIES_MAG2contigs = ();
open IN, "/storage1/data11/TYMEFLIES_phage/Binning_Data/TYMEFLIES_all_MAGs_stat.txt";
while (<IN>){
	chomp;
	
}
close IN;

## Step 4.2 Store all contigs that are determined to be viral sequence mis-binning
          # And also make filtered fsa files and concatenate them into a fasta and make blastn db based on it
my %Contig_viral_seq_in_MAGs = (); # $contig => $mag (the MAG that contains the contig)
open IN, "ls /storage1/data11/TYMEFLIES_phage/Binning_Data/split_fsa/*.fsa |";
while (<IN>){
	chomp;
	my $fsa_file = $_;
	my ($out_name) = $fsa_file =~ /split_fsa\/(.+?)\.fsa/;
	my %Seq = _store_seq("$fsa_file");
	my %Contig2len = (); # $contig => $len
	foreach my $key (sort keys %Seq){
		my ($key_clean) = $key =~ /^>(.+?)$/;
		my $len = length($Seq{$key});
		$Contig2len{$key_clean} = $len;
	}
	
	# Read the blastn out
	open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/filter_MAGs_out/TYMEFLIES_MAGs_${out_name}.blastn_out.txt";
	while (<IN>){
		chomp;
		my @tmp = split (/\t/);
		my $contig = $tmp[0];
		my $qstart = $tmp[6];
		my $qend = $tmp[7];
		my $viral_covered_len = abs($qend - $qstart);
		my $viral_covered_prec = $viral_part_len / $Contig2len{$contig};
		if ($viral_covered_prec >= 0.5){ # Cutoff: viral sequence hit covered region ≥ 50% of the host contig
			$contig
		}
	}
	close IN;
}

open IN, "ls /slowdata/databases/GEM/split_fsa/*.fsa |";
while (<IN>){
}





# Subroutine

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







