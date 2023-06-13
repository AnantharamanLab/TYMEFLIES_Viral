#!/usr/bin/perl

use strict;
use warnings;

# AIM: Run blastn to filter viral sequences from MAGs (both TYMEFLIES and GEM MAGs)
=pod
# Step 1. concatenate all nt sequences from both databases and split them into 20000 seq-containing fsa files
my $work_dir_TYMEFLIES_MAGs = "/storage1/data11/TYMEFLIES_phage/Robin_MAGs";
`mkdir $work_dir_TYMEFLIES_MAGs/tmp_fasta_files`;
## Only concatenate TYMEFLIES MAGs with >= 50% completeness and < 10% contamination
my %TYMEFLIES_MAGs_all_seq = (); # Store all TYMEFLIES MAGs with >= 50% completeness and < 10% contamination
open IN, "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.txt";
while (<IN>){
	chomp;
	if (!/^tymeflies/){
		my @tmp = split (/\t/);
		my $bin_full_name = $tmp[5];
		my $num_in_cluster = $tmp[15];
		if ($num_in_cluster ne "NA"){
			my ($img) = $bin_full_name =~ /_(33\d+?)_/;
			my $MAG_addr = $work_dir_TYMEFLIES_MAGs."/".$img."/".$bin_full_name.".fasta";
			`cp $MAG_addr $work_dir_TYMEFLIES_MAGs/tmp_fasta_files`;
		}
	}
}
close IN;

`find  $work_dir_TYMEFLIES_MAGs/tmp_fasta_files -name '*.fasta' -exec cat {} + > $work_dir_TYMEFLIES_MAGs/TYMEFLIES_MAGs_all_seq.fna`;
`rm -rf $work_dir_TYMEFLIES_MAGs/tmp_fasta_files`;

my $work_dir_GEM_MAGs = "/storage1/databases/GEM";
`find  $work_dir_GEM_MAGs/genomes -name '*.fasta' -exec cat {} + > $work_dir_GEM_MAGs/GEM_MAGs_all_seq.fna`;

`mkdir $work_dir_TYMEFLIES_MAGs/split_fsa`;
`perl /storage1/data11/TYMEFLIES_phage/split_multifasta.pl --in $work_dir_TYMEFLIES_MAGs/TYMEFLIES_MAGs_all_seq.fna --output_dir=$work_dir_TYMEFLIES_MAGs/split_fsa --seqs_per_file=20000`;

`mkdir $work_dir_GEM_MAGs/split_fsa`;
`perl /storage1/data11/TYMEFLIES_phage/split_multifasta.pl --in $work_dir_GEM_MAGs/GEM_MAGs_all_seq.fna --output_dir=$work_dir_GEM_MAGs/split_fsa --seqs_per_file=20000`;

`rm $work_dir_TYMEFLIES_MAGs/TYMEFLIES_MAGs_all_seq.fna $work_dir_GEM_MAGs/GEM_MAGs_all_seq.fna`;


# Step 2. Make viral sequence databases (blastn database)
# Viral sequence databases include: 1) IMG VR v4 all phages 2) NCBI RefSeq all viruses
#                                   3) TYMEFLIES all phages 4) Lake Mendota time series (2008-2012) all phages 

`mkdir All_viral_seq_db`;

#/storage1/databases/IMGVR-NCBI_phages/IMGVR_V4/IMGVR_all_nucleotides.fna  -> 15,722,824 sequences
#/storage1/databases/NCBI_RefSeq_viral/viral.genomic.2023-03-13.fna -> 15,288 sequences
`find /storage1/data11/TYMEFLIES_phage/33*/vRhyme_best_bins_fasta_parsed/ -name '*.fasta' -exec cat {} + > All_viral_seq_db/TYMEFLIES_all_phages.fasta`;  # 1,820,639 sequences
`find /storage1/data11/LakeMendota_2008_to_2012/33*/vRhyme_best_bins_fasta_parsed/ -name '*.fasta' -exec cat {} + > All_viral_seq_db/Lake_Mendota_2008_to_2012_all_phages.fasta`; # 34,071 sequences

`cat /storage1/databases/IMGVR-NCBI_phages/IMGVR_V4/IMGVR_all_nucleotides.fna /storage1/databases/NCBI_RefSeq_viral/viral.genomic.2023-03-13.fna All_viral_seq_db/TYMEFLIES_all_phages.fasta All_viral_seq_db/Lake_Mendota_2008_to_2012_all_phages.fasta > All_viral_seq_db/All_viral_seq.fasta`;
# 17,592,822 sequences in All_viral_seq_db/All_viral_seq.fasta

`rm All_viral_seq_db/Lake_Mendota_2008_to_2012_all_phages.fasta`;

`makeblastdb -in All_viral_seq_db/All_viral_seq.fasta -title All_viral_seq -dbtype nucl -out All_viral_seq_db/All_viral_seq_blastdb`;

# Step 3. Run blastn 
# Step 3.1 Write down blastn command for TYMEFLIES MAGs and GEM MAGs
`mkdir Host_prediction`;
`mv All_viral_seq_db/TYMEFLIES_all_phages.fasta Host_prediction/All_phage_genomes.fasta`;
`mkdir Host_prediction/filter_MAGs_out`;

open OUT, ">tmp.filter_MAGs_blastn_1.sh";
open IN, "ls /storage1/data11/TYMEFLIES_phage/Robin_MAGs/split_fsa/*.fsa |";
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
open IN, "ls /storage1/databases/GEM/split_fsa/*.fsa |";
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

# Step 3.2 Run the command with 20 cpus

`cat tmp.filter_MAGs_blastn_1.sh | parallel -j 20`;
`cat tmp.filter_MAGs_blastn_2.sh | parallel -j 20`;

`rm tmp.filter_MAGs_blastn_1.sh`;
`rm tmp.filter_MAGs_blastn_2.sh`;
=cut
my $work_dir_TYMEFLIES_MAGs = "/storage1/data11/TYMEFLIES_phage/Robin_MAGs";
# Step 4. Parse result to get viral genome to final host prediction
## Step 4.1 Store all contig to MAG info and MAG to tax info
my %TYMEFLIES_contig2MAG = (); # $contig => $mag
my %TYMEFLIES_MAG2contigs = (); # $mag => $contig collection separated by "\,"
my %MAG2tax = (); # $mag => $tax
open IN, "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.txt";
while (<IN>){
	chomp;
	if (!/^tymeflies/){
		my @tmp = split (/\t/);
		my $mag = $tmp[5];
		my $num_in_cluster = $tmp[15];		
		if ($num_in_cluster ne "NA"){
			my ($img) = $mag =~ /_(33\d+?)_/;
			my @Contigs = (); # Store all the contigs into an array
			my $MAG_addr = $work_dir_TYMEFLIES_MAGs."/".$img."/".$mag.".fasta";
			my %MAG_seq = _store_seq("$MAG_addr");
			foreach my $header (sort keys %MAG_seq){
				my ($contig) = $header =~ /^>(.+?)$/;
				push @Contigs, $contig;
				$TYMEFLIES_contig2MAG{$contig} = $mag;
			}
			
			$TYMEFLIES_MAG2contigs{$mag} = join(",", @Contigs);
			my $tax = join(";", @tmp[16..22]);
			$MAG2tax{$mag} = $tax;
		}
	}
}
close IN;

my %GEM_contig2MAG = (); # $contig => $mag
my %GEM_MAG2contigs = (); # $mag => $contig collection separated by "\,"
open IN, "/storage1/databases/GEM/GEM_all_MAGs_stat.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $contigs = $tmp[4];
		my $tax = $tmp[1];
		$GEM_MAG2contigs{$mag} = $contigs;
		
		my @Contigs = split (/\,/, $contigs);
		foreach my $contig (@Contigs){
			$GEM_contig2MAG{$contig} = $mag;
		}
		
		$MAG2tax{$mag} = $tax;		
	}
}
close IN;

my %Contig2MAG = (%TYMEFLIES_contig2MAG, %GEM_contig2MAG);
my %MAG2contigs = (%TYMEFLIES_MAG2contigs, %GEM_MAG2contigs);

## Step 4.2 Store all contigs that are determined to be viral sequence mis-binning,
#           and store the match result of viral sequence to contig   
my %Contig_viral_seq_in_MAGs = (); # $contig => $mag (the MAG that contains the contig)
my %Viral_seq2contig = (); # $viral_seq => $contig collection separated by "\,"
open IN, "ls /storage1/data11/TYMEFLIES_phage/Robin_MAGs/split_fsa/*.fsa |";
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
	open INN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/filter_MAGs_out/TYMEFLIES_MAGs_${out_name}.blastn_out.txt";
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $contig = $tmp[0];
		my $viral_seq = $tmp[1];
		my $qstart = $tmp[6];
		my $qend = $tmp[7];
		my $tstart = $tmp[8];		
		my $tend = $tmp[9];		
		my $iden = $tmp[2];
		my $host_covered_len = abs($qend - $qstart);
		my $viral_covered_len = abs($tend - $tstart);		
		my $host_covered_prec = $host_covered_len / $Contig2len{$contig};
		if ($host_covered_prec >= 0.5 and $iden >= 80){ # Cutoff: viral sequence hit covered region ≥ 50% of the host contig
			my $mag = $Contig2MAG{$contig};
			$Contig_viral_seq_in_MAGs{$contig} = $mag;
		}
		if ($host_covered_len >= 2000 and $viral_covered_len >= 2000 and $iden >= 95){ #  Host predictions were then based on matches of ≥90% nucleotide identity covering ≥2 kb of the virus and (putative) host sequences.
			if (! exists $Viral_seq2contig{$viral_seq}){
				$Viral_seq2contig{$viral_seq} = $contig;
			}else{
				$Viral_seq2contig{$viral_seq} .= "\,".$contig;
			}
		}
	}
	close INN;
}
close IN;

open IN, "ls /storage1/databases/GEM/split_fsa/*.fsa |";
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
	open INN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/filter_MAGs_out/GEM_MAGs_${out_name}.blastn_out.txt";
	while (<INN>){
		chomp;
		my @tmp = split (/\t/);
		my $contig = $tmp[0];
		my $viral_seq = $tmp[1];
		my $qstart = $tmp[6];
		my $qend = $tmp[7];
		my $tstart = $tmp[8];		
		my $tend = $tmp[9];		
		my $iden = $tmp[2];
		my $host_covered_len = abs($qend - $qstart);
		my $viral_covered_len = abs($tend - $tstart);		
		my $host_covered_prec = $host_covered_len / $Contig2len{$contig};
		if ($host_covered_prec >= 0.5 and $iden >= 80){ # Cutoff: viral sequence hit covered region ≥ 50% of the host contig
			my $mag = $Contig2MAG{$contig};
			$Contig_viral_seq_in_MAGs{$contig} = $mag;
		}
		if ($host_covered_len >= 2000 and $viral_covered_len >= 2000 and $iden >= 95){ #  Host predictions were then based on matches of ≥90% nucleotide identity covering ≥2 kb of the virus and (putative) host sequences.
			if (! exists $Viral_seq2contig{$viral_seq}){
				$Viral_seq2contig{$viral_seq} = $contig;
			}else{
				$Viral_seq2contig{$viral_seq} .= "\,".$contig;
			}
		}
	}
	close INN;
}
close IN;

## Step 4.3 Store the match result of viral sequence to host genome tax  
my %Viral_seq2host_tax = (); # $viral_seq => $host_tax collection separated by "\t"
foreach my $viral_seq (sort keys %Viral_seq2contig){
	my @Contigs = split (/\,/, $Viral_seq2contig{$viral_seq});
	my @Host_tax = (); # Host tax collection
	
	foreach my $contig (@Contigs){
		if (! exists $Contig_viral_seq_in_MAGs{$contig}){
			my $mag = $Contig2MAG{$contig};
			if (exists $TYMEFLIES_MAG2contigs{$mag}){ # Only use the genome from TYMEFLIES
				my $tax = $MAG2tax{$mag};
				push @Host_tax, $tax;
			}
		}
	}
	
	if (@Host_tax){ 
		$Viral_seq2host_tax{$viral_seq} = join("\t", @Host_tax);
	}
}

# TEST
# Print %Viral_seq2host_tax
open OUT, ">Host_prediction/Viral_seq2host_tax.txt";
foreach my $viral_seq (sort keys %Viral_seq2host_tax){
	print OUT "$viral_seq\t$Viral_seq2host_tax{$viral_seq}\n";
}
close OUT;

## Step 4.4 Store the match result of viral genome to host genome tax 
### Store all viral genome to viral sequences hash
my %Viral_gn2viral_seq = (); # $viral_gn => $viral_seq collection separated by "\t"
open IN, "/storage1/data11/TYMEFLIES_phage/Host_prediction/All_phage_genomes.fasta";
while (<IN>){
	chomp;
	if (/^>/){
		my $line = $_;
		my ($viral_seq) = $line =~ /^>(.+?)$/;
		my ($viral_gn) = $viral_seq =~ /^(.+?\_\_.+?)\_\_.+?$/;
		if (! exists $Viral_gn2viral_seq{$viral_gn}){
			$Viral_gn2viral_seq{$viral_gn} = $viral_seq;
		}else{
			$Viral_gn2viral_seq{$viral_gn} .= "\t".$viral_seq;
		}
	}
}
close IN;

my %Viral_gn2host_tax = (); # $viral_gn => $host_tax collection separated by "\t"
foreach my $viral_gn (sort keys %Viral_gn2viral_seq){
	my @Host_tax = ();
	my @Viral_seqs = split (/\t/, $Viral_gn2viral_seq{$viral_gn});
	foreach my $viral_seq (@Viral_seqs){
		if (exists $Viral_seq2host_tax{$viral_seq}){
			push @Host_tax, $Viral_seq2host_tax{$viral_seq};
		}
	}
	
	if (@Host_tax){
		$Viral_gn2host_tax{$viral_gn} = join("\t", @Host_tax);
	}
}

## Step 4.5 Get the final host ranks
my %Viral_gn2final_host_tax = (); # $viral_gn => $final_host_tax (based on 80% consensus rule at each rank)
foreach my $viral_gn (sort keys %Viral_gn2host_tax){
	my @Host_tax = split (/\t/, $Viral_gn2host_tax{$viral_gn});
	my $final_host_tax = _find_consensus_lineages_based_on_each_rank(@Host_tax);
	if ($final_host_tax){
		$Viral_gn2final_host_tax{$viral_gn} = $final_host_tax;
	}
}

# Step 4.6 Write down result
open OUT, ">Host_prediction/Phage_gn2host_tax_by_sequence_similarity.txt";
foreach my $gn (sort keys %Viral_gn2final_host_tax){
	print OUT "$gn\t$Viral_gn2final_host_tax{$gn}\n";
}
close OUT;



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

sub _find_consensus_lineages_based_on_each_rank{ # The default is 80% consensus
	my @Lineages = @_; # Get the passed array
	
	my %Domain = (); # $domain => the times of this domain appears
	my %Phylum = (); # $phylum => the times of this phylum appears
	my %Class = (); # $class => the times of this class appears
	my %Order = (); # $order => the times of this order appears
	my %Family = (); # $family => the times of this family appears
	my %Genus = (); # $genus => the times of this genus appears
	my %Species = (); # $species => the times of this species appears
	my $lineage_hit_num = scalar @Lineages;
	
	foreach my $lineage (@Lineages){
		my ($domain,$phylum,$class,$order,$family,$genus,$species) = $lineage =~ /^(.+?)\;(.+?)\;(.+?)\;(.+?)\;(.+?)\;(.+?)\;(.+?)$/;
		$Domain{$domain}++;
		$Phylum{$phylum}++;
		$Class{$class}++;
		$Order{$order}++;
		$Family{$family}++;
		$Genus{$genus}++;
		$Species{$species}++;
	}
	
	my $lineage_final = "";
	my @Lineage_final = ();
	
	foreach my $domain (sort keys %Domain){
		if (($Domain{$domain} / $lineage_hit_num) >= 0.8){
			$Lineage_final[0] = $domain;
		}else{
			$Lineage_final[0] = "Not consensual";
		}
	}
		
	foreach my $phylum (sort keys %Phylum){
		if (($Phylum{$phylum} / $lineage_hit_num) >= 0.8){
			$Lineage_final[1] = $phylum;
		}else{
			$Lineage_final[1] = "Not consensual";
		}
	}

	foreach my $class (sort keys %Class){
		if (($Class{$class} / $lineage_hit_num) >= 0.8){
			$Lineage_final[2] = $class;
		}else{
			$Lineage_final[2] = "Not consensual";
		}
	}
	
	foreach my $order (sort keys %Order){
		if (($Order{$order} / $lineage_hit_num) >= 0.8){
			$Lineage_final[3] = $order;
		}else{
			$Lineage_final[3] = "Not consensual";
		}
	}	
	
	foreach my $family (sort keys %Family){
		if (($Family{$family} / $lineage_hit_num) >= 0.8){
			$Lineage_final[4] = $family;
		}else{
			$Lineage_final[4] = "Not consensual";
		}
	}	

	foreach my $genus (sort keys %Genus){
		if (($Genus{$genus} / $lineage_hit_num) >= 0.8){
			$Lineage_final[5] = $genus;
		}else{
			$Lineage_final[5] = "Not consensual";
		}
	}	

	foreach my $species (sort keys %Species){
		if (($Species{$species} / $lineage_hit_num) >= 0.8){
			$Lineage_final[6] = $species;
		}else{
			$Lineage_final[6] = "Not consensual";
		}
	}	
	
	my @Lineage_final_curated = (); # Curated final lineage
	for(my $i=0; $i<=$#Lineage_final; $i++){
		if ($Lineage_final[$i] ne "Not consensual" and $Lineage_final[$i] !~ /^\S\_\_$/){
			push @Lineage_final_curated, $Lineage_final[$i];
		}else{
			last;
		}
	}
	
	my $lineage_final_curated = join("\;",@Lineage_final_curated);
	
	return $lineage_final_curated;
}





