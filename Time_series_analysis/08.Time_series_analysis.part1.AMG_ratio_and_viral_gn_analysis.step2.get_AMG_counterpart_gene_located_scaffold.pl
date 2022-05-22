#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get AMG counterpart gene located scaffolds from all metagenome scaffolds

# Step 1 Store all AMG info
my %AMG_pro2ko = (); # $amg_pro => $ko
my %KO2amg_pro = (); # $ko => $amg_pros (collection of $amg_pro, separated by "\,")
open IN, "/storage1/data11/TYMEFLIES_phage/AMG_analysis/AMG_summary.txt";
while (<IN>){
	chomp;
	if (!/^Pro/){
		my @tmp = split (/\t/);
		my $amg_pro = $tmp[0];
		my $ko = $tmp[2];
		$AMG_pro2ko{$amg_pro} = $ko;
		
		if (!exists $KO2amg_pro{$ko}){
			$KO2amg_pro{$ko} = $amg_pro;
		}else{
			$KO2amg_pro{$ko} .= "\,".$amg_pro;
		}
	}
}
close IN;

# Step 2 Find AMG counterpart gene located scaffolds
## Step 2.1 Perform hmmsearch for all KO in %KO2amg_pro on all metagenome proteins
### Step 2.1.1 Store KO2cutoff_and_type hash
my %KO2cutoff_and_type = (); # $ko => $cutoff_and_type
open IN, "/slowdata/data1/Genome_profile_software/kofam_database/ko_list";
while (<IN>){
	chomp;
	if (!/^knum/){
		my @tmp = split (/\t/);
		my $ko = $tmp[0];
		my $cutoff_and_type = $tmp[1]."\t".$tmp[2];
		$KO2cutoff_and_type{$ko} = $cutoff_and_type;
	}
}
close IN;

### Step 2.1.2 Store all IMG ID
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '^33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;
=pod
### Step 2.1.3 Run hmmsearch
`mkdir tmp_folder_for_hmmsearch_for_AMG_counterpart_genes`;
open OUT, ">tmp.run_hmmsearch_for_getting_AMG_counterpart_gene_located_scaffold.sh";
foreach my $img_id (sort keys %IMGID){
	foreach my $ko (sort keys %KO2amg_pro){
		if (exists $KO2cutoff_and_type{$ko}){
			my ($cutoff, $type) = $KO2cutoff_and_type{$ko} =~ /^(.+?)\t(.+?)$/;
			if ($cutoff eq "\-"){
				$cutoff = 50;
			}
			if (!-e ("tmp_folder_for_hmmsearch_for_AMG_counterpart_genes/$img_id.$ko.hmmsearch_result.txt")){
				if ($type eq "full"){
					print OUT "hmmsearch -T $cutoff --cpu 1 --tblout tmp_folder_for_hmmsearch_for_AMG_counterpart_genes/$img_id.$ko.hmmsearch_result.txt /slowdata/data1/Genome_profile_software/kofam_database/profiles/$ko.hmm $img_id/$img_id.a.faa\n";
				}else{
					print OUT "hmmsearch --domT $cutoff --cpu 1 --tblout tmp_folder_for_hmmsearch_for_AMG_counterpart_genes/$img_id.$ko.hmmsearch_result.txt /slowdata/data1/Genome_profile_software/kofam_database/profiles/$ko.hmm $img_id/$img_id.a.faa\n";
				}
			}
		}
	}
}

`cat tmp.run_hmmsearch_for_getting_AMG_counterpart_gene_located_scaffold.sh | parallel -j 30`;

`rm tmp.run_hmmsearch_for_getting_AMG_counterpart_gene_located_scaffold.sh`;
=cut
## Step 2.2 Read and parse hmmsearch results
my %IMG2all_pros = (); # $img_id => $all_pros (collection of $pro, separated by "\,"); all protein hits of AMG KOs
open IN, "find /storage1/data11/TYMEFLIES_phage/tmp_folder_for_hmmsearch_for_AMG_counterpart_genes/ -name '*.hmmsearch_result.txt' | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($img_id, $ko) = $file =~ /^.+\/(33\d+?)\.(K.+?)\.hmmsearch_result\.txt/;
	my ($cutoff, $type) = $KO2cutoff_and_type{$ko} =~ /^(.+?)\t(.+?)$/;

	if ($cutoff eq "\-"){
		$cutoff = 50;
	}	
	
	open INN, "$file";
	while (<INN>){
		chomp;
		if (!/^#/){
			my $line = $_;
			$line =~ s/ +/ /g;
			my @tmp = split (/\s/, $line);
			my $pro = $tmp[0];
			my $score_full = $tmp[5];
			my $score_domain = $tmp[8];
			if ($type eq "full"){
				if ($score_full >= $cutoff){
					if (!exists $IMG2all_pros{$img_id}){
						$IMG2all_pros{$img_id} = $pro;
					}else{
						$IMG2all_pros{$img_id} .= "\,".$pro;
					}
				}
			}else{
				if ($score_domain >= $cutoff){
					if (!exists $IMG2all_pros{$img_id}){
						$IMG2all_pros{$img_id} = $pro;
					}else{
						$IMG2all_pros{$img_id} .= "\,".$pro;
					}				
				}
			}
		}
	}
	close INN;
}
close IN;

## Step 2.3 Get AMG counterpart gene located scaffolds
my %Scaffolds_from_viruses = (); # $scaffold => 1
open IN, "Host_prediction/All_phage_genomes_headers.txt";
while (<IN>){
	chomp;
	my ($scaffold) = $_ =~ /^>.+?\_\_(Ga.+?)$/;
	$Scaffolds_from_viruses{$scaffold} = 1;
}
close IN;

my %AMG_counterpart_gene_located_scaffolds = (); # $img_id => $scaffolds (collection of $scaffold, separated by "\,")
foreach my $img_id (sort keys %IMG2all_pros){
	# Store all scaffolds ID that contain AMG KO hits
	my %Scaffolds_all = (); # $scaffold => 1
	my @Pros = split (/\,/, $IMG2all_pros{$img_id});
	foreach my $pro (@Pros){
		my ($scaffold) = $pro =~ /^(Ga.+?\_.+?)\_/;
		if (! exists $Scaffolds_from_viruses{$scaffold}){ # Not include those scaffolds from viral genomes
			$Scaffolds_all{$scaffold} = 1;
		}
	}
	
	my @Scaffolds_all = sort keys (%Scaffolds_all);
	my $scaffolds = join("\,", @Scaffolds_all);
	
	$AMG_counterpart_gene_located_scaffolds{$img_id} = $scaffolds; 
}

## Step 2.4 Write down AMG_counterpart_gene_located_scaffolds.fasta
my %Seq_for_AMG_counterpart_gene_located_scaffolds = (); 
foreach my $img_id (sort keys %AMG_counterpart_gene_located_scaffolds){
	# Store all the scaffold sequences
	my %Seq = _store_seq("$img_id/$img_id.a.fna");
	
	# Add AMG counterpart gene located scaffolds into the hash
	my @Scaffolds = split (/\,/, $AMG_counterpart_gene_located_scaffolds{$img_id});
	
	foreach my $scaffold (@Scaffolds){
		my $scaffold_w_array = ">".$scaffold;
		$Seq_for_AMG_counterpart_gene_located_scaffolds{$scaffold_w_array} = $Seq{$scaffold_w_array};
	}
}

open OUT, ">AMG_counterpart_gene_located_scaffolds.fasta";
foreach my $key (sort keys %Seq_for_AMG_counterpart_gene_located_scaffolds){
	print OUT "$key\n$Seq_for_AMG_counterpart_gene_located_scaffolds{$key}\n";
}
close OUT;



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