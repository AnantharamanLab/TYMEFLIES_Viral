#!/usr/bin/perl

use strict;
use warnings;
use Time::HiRes qw(usleep nanosleep);

# Aim: Cluster phage genomes by share protein and AAI
=pod
# Step 1 Make diamond db to all phage genome proteins
## Concatenate all phage genome proteins
`mkdir Cluster_phage_genomes`;
`find /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed -name '*.faa' | xargs cat > Cluster_phage_genomes/All_phage_genome.faa`;
# All_phage_genome.faa has 18,578,326 sequences

## Make db
`diamond makedb --in Cluster_phage_genomes/All_phage_genome.faa --db Cluster_phage_genomes/All_phage_genome_proteins --threads 15`;

# Step 2 Perform all-vs-all BLASTP
`mkdir Cluster_phage_genomes/All_phage_genome_faa_split_fsa`;
`perl /storage1/data11/TYMEFLIES_phage/split_multifasta.pl --in Cluster_phage_genomes/All_phage_genome.faa --output_dir=Cluster_phage_genomes/All_phage_genome_faa_split_fsa --seqs_per_file=40000`;

open OUT, ">tmp.run_all_vs_all_diamond_blastp.sh";
open IN, "ls Cluster_phage_genomes/All_phage_genome_faa_split_fsa/*.fsa |";
while (<IN>){
	chomp;
	my $fsa_full_adr = $_;
	my ($fsa_name) = $fsa_full_adr =~ /All_phage_genome_faa_split_fsa\/(.+?)\.fsa/;
	print OUT "diamond blastp --query $fsa_full_adr --db Cluster_phage_genomes/All_phage_genome_proteins --out Cluster_phage_genomes/All_phage_genome.$fsa_name.diamond_blastp.tsv --outfmt 6 --evalue 1e-5 --max-target-seqs 10000 --query-cover 50 --subject-cover 50 --threads 1\n";

}
close IN;
close OUT;

`cat tmp.run_all_vs_all_diamond_blastp.sh | parallel -j 15`;
`rm tmp.run_all_vs_all_diamond_blastp.sh`;
=cut
# Step 3 Parse the result of all-vs-all diamond blastp
## Step 3.1 Store bin to bin size (protein number)
my %Bin2bin_size = (); # $bin => $bin_size
open IN, "find /storage1/data11/TYMEFLIES_phage/*/vRhyme_best_bins_fasta_parsed/Each_bin_info.txt |";
while (<IN>){
	chomp;
	my $file = $_;
	open INN, $file;
	while (<INN>){
		chomp;
		if (!/^bin/){
			my @tmp = split (/\t/);
			my $bin = $tmp[0];
			my $bin_size = $tmp[3];
			$Bin2bin_size{$bin} = $bin_size;
		}
	}
	close INN;
}
close IN;

## Step 3.2 Store protein to bin map
my %Pro2bin = (); # $pro => $bin
open IN, "find 33* -maxdepth 1 -type d -name 'vRhyme_best_bins_fasta_parsed' |";
while (<IN>){
	chomp;
	my $folder_adr = $_;
	open INN, "grep '>' $folder_adr/*.faa |";
	while (<INN>){
		chomp;
		my ($bin,$pro) = $_ =~ /vRhyme_best_bins_fasta_parsed\/(.+?)\.faa\:\>(.+?)$/;
		$Pro2bin{$pro} = $bin;
	}
	close INN;
}
close IN;

## Step 3.3 Make bin vs bin hash for comparison
#`cat Cluster_phage_genomes/All_phage_genome.*.diamond_blastp.tsv > Cluster_phage_genomes/All_phage_genome_diamond_blastp.tsv`;
#`rm Cluster_phage_genomes/All_phage_genome.*.diamond_blastp.tsv`;
=pod
my %Bin_list = (); # $bin => $bin_order
my $i = 1;
foreach my $bin (sort keys %Bin2bin_size){
	$Bin_list{$bin} = $i;
	$i++;
}

`mkdir Cluster_phage_genomes/Bin_vs_bin_comparison`;

open IN, "ls Cluster_phage_genomes/All_phage_genome.*.diamond_blastp.tsv |";
while (<IN>){
	chomp;
	my $diamond_blastp_result = $_;
	
	my %Hash = (); # For each All_phage_genome.*.diamond_blastp.tsv, store the result name and result content
				   # $file => $aln collection separated by "\t"
	open INN, "$diamond_blastp_result";
	while (<INN>){
		chomp;
		my $aln = $_;
		my @tmp = split (/\t/,$aln);
		my $query = $Pro2bin{$tmp[0]};
		my $target = $Pro2bin{$tmp[1]};
		my $query_order = $Bin_list{$query};
		my $target_order = $Bin_list{$target};
		
		if ($query ne $target and $query_order < $target_order){
			my $file = "Cluster_phage_genomes/Bin_vs_bin_comparison/$query\_\_vs\_\_others.diamond_blastp.tsv";
			if (! exists $Hash{$file}){
				$Hash{$file} = $aln;
			}else{
				$Hash{$file} .= "\,".$aln;
			}
		}
	}
	close INN;
	
	foreach my $file (sort keys %Hash){
		my @Aln = split (/\,/, $Hash{$file});
		open OUT, ">>$file";
		# 1 millisecond == 1000 microseconds
		usleep(10000); # Sleep 10,000 microseconds, which is 10 millisecond, 0.01 second	
		foreach my $aln (@Aln){
			print OUT "$aln\n";
		}
		close OUT;
	}
}
close IN;
=cut
## Step 3.4 Calculate aai for each file
`mkdir /fastdata/archive/Shared_protein_and_AAI_result`;
open IN, "find Cluster_phage_genomes/Bin_vs_bin_comparison -name '*.diamond_blastp.tsv' |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($file_name) = $file =~ /Bin_vs_bin_comparison\/(.+?)\.diamond_blastp\.tsv/;
	my %Hash2aai_n_score = (); # Store the alignment result for each pair of comparison, the best hit of each gene
	                           # $query => $target => $gene => $aai_n_score
	my %Query = (); # $query => 1
	my %Target = (); # $target => 1
	open INN, $file;
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $gene = $tmp[0];
		my $aai = $tmp[2];
		my $score = $tmp[-1];
		my $query = $Pro2bin{$tmp[0]};
		my $target = $Pro2bin{$tmp[1]};
		
		$Query{$query} = 1;
		$Target{$target} = 1;
		
		if (! exists $Hash2aai_n_score{$query}{$target}{$gene}){
			$Hash2aai_n_score{$query}{$target}{$gene} = "$aai\t$score";
		}else{
			my ($aai_previous, $score_previous) = $Hash2aai_n_score{$query}{$target}{$gene} =~ /^(.+?)\t(.+?)$/;
			if ($score > $score_previous){
				$Hash2aai_n_score{$query}{$target}{$gene} = "$aai\t$score";
			}
		}
	}
	close INN;
	
	my @Row = (); # Store the result
	foreach my $query (sort keys %Query){
		foreach my $target (sort keys %Target){
			my $query_genes = $Bin2bin_size{$query}; 
			my $target_genes = $Bin2bin_size{$target};
			
			my $shared_genes = 0;
			my $aai = 0; # The average AAI for all shared genes
			my $aai_sum = 0; # The sum AAI for all shared genes
			
			foreach my $gene (sort keys %{$Hash2aai_n_score{$query}{$target}}){
				$shared_genes++;
				my ($aai,$score) = $Hash2aai_n_score{$query}{$target}{$gene} =~ /^(.+?)\t(.+?)$/;
				$aai_sum += $aai;
			}
			
			my $qcov = 0; 
			my $tcov = 0;
			if ($shared_genes){
				$aai = $aai_sum / $shared_genes; 
				$qcov = $shared_genes / $query_genes * 100;
				$tcov = $shared_genes / $target_genes * 100;
			}	
			
			my $row = "$query\t$target\t$query_genes\t$target_genes\t$shared_genes\t$qcov\t$tcov\t$aai";	
			push @Row, $row;
		}
	}
	
	# Write down the result
	open OUT, ">/fastdata/archive/Shared_protein_and_AAI_result/Shared_protein_and_AAI.$file_name.txt";
	foreach my $row (@Row){
		print OUT "$row\n";
	}
	close OUT;
	
}
close IN;
