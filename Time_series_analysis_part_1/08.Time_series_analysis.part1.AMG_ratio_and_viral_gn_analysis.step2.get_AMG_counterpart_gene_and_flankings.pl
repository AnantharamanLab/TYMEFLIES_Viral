#!/usr/bin/perl

use strict;
use warnings;

# Aim: Get AMG counterpart genes from all metagenomes
# Note there are two requirements: 
# 1) Only include hits from microbial scaffolds
# 2) Also include flanking regions (150 bp on left and right)

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

`cat tmp.run_hmmsearch_for_getting_AMG_counterpart_gene_located_scaffold.sh | parallel -j 20`;

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

open OUT, ">IMG2all_pros.txt";
foreach my $img_id (sort keys %IMG2all_pros){
	print OUT "$img_id\t$IMG2all_pros{$img_id}\n";
}
close OUT;

## Step 2.3 Get AMG counterpart gene and flanking regions
my %Scaffolds_from_viruses = (); # $scaffold => 1; do not contain "fragment" in the header
open IN, "Host_prediction/All_phage_genomes_headers.txt";
while (<IN>){
	chomp;
	my ($scaffold) = $_ =~ /^>.+?\_\_(Ga.+?\_.+?)$/;
	if ($scaffold =~ /fragment/){
		my ($scaffold_wo_fragment) = $scaffold =~ /^(.+)_fragment/;
		$Scaffolds_from_viruses{$scaffold_wo_fragment} = 1;
	}else{
		$Scaffolds_from_viruses{$scaffold} = 1;
	}
	
}
close IN;

my %Scaffold2length = (); # $scaffold => $length
open IN, "ls */*seq_length.txt |";
while (<IN>){
	chomp;
	my $file = $_;
	open INN, $file;
	while (<INN>){
		chomp;
		if (!/^sequence/){
			my @tmp = split (/\t/);
			my $scaffold = $tmp[0];
			my $length = $tmp[1];
			$Scaffold2length{$scaffold} = $length;
		}
	}
	close INN;
}
close IN;

my %Seq_for_AMG_counterpart_genes_and_flankings = (); # Store all the sequence 
my %All_scaffold2gene_and_flanking_regions = (); # $scaffold => $regions; for all scaffolds
foreach my $img_id (sort keys %IMG2all_pros){
	# Store all the scaffold sequences
	my %Seq = _store_seq("$img_id/$img_id.a.fna");

	# Get the scaffold to protein ID
	my @Pros = split (/\,/, $IMG2all_pros{$img_id});
	my %Scaffold2pros = (); # Store the proteins that are only from microbial scaffolds
	foreach my $pro (@Pros){
		my ($scaffold) = $pro =~ /^(Ga.+?\_.+?)\_/;
		if (! exists $Scaffolds_from_viruses{$scaffold}){ # Not include those scaffolds from viral genomes
			if (! exists $Scaffold2pros{$scaffold}){
				$Scaffold2pros{$scaffold} = $pro;
			}else{
				$Scaffold2pros{$scaffold} .= "\,".$pro;
			}
		}
	}	
	
	# Get scaffold to gene and flanking regions for each scaffold and store the sequences
	foreach my $scaffold (sort keys %Scaffold2pros){
		my $length = $Scaffold2length{$scaffold};
		
		my @Regions = (); # Store the gene and flanking regions
		my @Pros_from_this_scaffold = split (/\,/, $Scaffold2pros{$scaffold});
		foreach my $pro (@Pros_from_this_scaffold){
			my ($start, $end) = $pro =~ /^Ga.+?\_.+?\_(.+?)\_(.+?)$/;
			if ($start >= 151){
				$start = $start - 150;
			}else{
				$start = 1;
			}
			
			if (($end + 150) <= $length){
				$end = $end + 150;
			}else{
				$end = $length;
			}
			
			push @Regions, "$start\_$end";
		}
		
		my @Regions_merged = _merge_regions(@Regions);
		$All_scaffold2gene_and_flanking_regions{$scaffold} = join("\,", @Regions_merged);
		
		# Store sequences
		foreach my $key (@Regions_merged){
			my ($start, $end) = $key =~ /^(.+?)\_(.+?)$/;
			my $head = ">$scaffold\_$start\_$end";
			my $seq = substr($Seq{">$scaffold"}, ($start - 1), ($end - $start + 1));
			if (! $seq){
				print "$head is empty!\n";
			}else{
				$Seq_for_AMG_counterpart_genes_and_flankings{$head} = $seq;
			}
		}
	}
	
	print "$img_id has been processed!\n";
}

open OUT, ">AMG_counterpart_genes_and_flankings.fasta";
foreach my $key (sort keys %Seq_for_AMG_counterpart_genes_and_flankings){
	print OUT "$key\n$Seq_for_AMG_counterpart_genes_and_flankings{$key}\n";
}
close OUT;

open OUT, ">All_scaffold2gene_and_flanking_regions.txt";
foreach my $key (sort keys %All_scaffold2gene_and_flanking_regions){
	print OUT "$key\t$All_scaffold2gene_and_flanking_regions{$key}\n";
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

sub _merge_regions{
	my @List = @_;
	my @Result = (); 
	
	my @Pair_data = (); # $start\_-1 or $end\_1
	for(my $i=0; $i<=$#List; $i++){
		my ($start, $end) = $List[$i] =~ /^(.+?)\_(.+?)$/;
		my $ele_1 = $start."\_-1";
		my $ele_2 = $end."\_1";
		push @Pair_data, $ele_1;
		push @Pair_data, $ele_2;
	}	
	
	# Sort the pair data in ascending order according to the first part of the element
	@Pair_data = sort { ($a =~ /^(\d+)_/)[0] <=> ($b =~ /^(\d+)_/)[0] } @Pair_data;
	
	my $indicator = 0;
	my $mark_start = -1;
	my $mark_end = -1;
	my $start_flagged = 0;
	
	for(my $i=0; $i<=$#Pair_data; $i++){
		$indicator += (split(/\_/, $Pair_data[$i]))[1]; 
		
		if ($indicator == -1 and $start_flagged == 0){ # The first entry point
			$mark_start = (split(/\_/, $Pair_data[$i]))[0];
			$start_flagged = 1;
		}elsif($indicator == 0){
			$mark_end = (split(/\_/, $Pair_data[$i]))[0];
			$start_flagged = 0;
		}
		
		if ($mark_end >= $mark_start){ # A pair has been found
			push @Result, "$mark_start\_$mark_end";
		}
	}
	
	return @Result;
}
