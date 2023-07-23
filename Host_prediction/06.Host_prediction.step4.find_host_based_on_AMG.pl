#!/usr/bin/perl

use strict;
use warnings;

# AIM: Find viral host based on the AMG identity 

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

# Step 2 KO to microbial pro hash
## Step 2.1 Perform hmmsearch for all KO in %KO2amg_pro
#`mkdir Host_prediction/Find_host_based_on_AMG`;
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
=pod
### Step 2.1.2 Run hmmsearch
my $all_seq = "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/TYMEFLIES_MAGs_all_seq.faa";
open OUT, ">tmp.run_hmmsearch_for_host_prediction.sh";
foreach my $ko (sort keys %KO2amg_pro){
	if (exists $KO2cutoff_and_type{$ko}){
		my ($cutoff, $type) = $KO2cutoff_and_type{$ko} =~ /^(.+?)\t(.+?)$/;
		if ($cutoff eq "\-"){
			$cutoff = 50;
		}
		if (!-e ("Host_prediction/Find_host_based_on_AMG/$ko.hmmsearch_result.txt")){
			if ($type eq "full"){
				print OUT "hmmsearch -T $cutoff --cpu 1 --tblout Host_prediction/Find_host_based_on_AMG/$ko.hmmsearch_result.txt /slowdata/data1/Genome_profile_software/kofam_database/profiles/$ko.hmm $all_seq\n";
			}else{
				print OUT "hmmsearch --domT $cutoff --cpu 1 --tblout Host_prediction/Find_host_based_on_AMG/$ko.hmmsearch_result.txt /slowdata/data1/Genome_profile_software/kofam_database/profiles/$ko.hmm $all_seq\n";
			}
		}
	}
}

`cat tmp.run_hmmsearch_for_host_prediction.sh | parallel -j 20`;

`rm tmp.run_hmmsearch_for_host_prediction.sh`;


### Step 2.1.3 Store all contig to MAG info and MAG to tax info
my %TYMEFLIES_contig2MAG = (); # $contig => $mag
open IN, "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.CheckM_passed.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $contigs = $tmp[4];
		my $tax = $tmp[1];
		my @Contigs = split (/\,/, $contigs);
		foreach my $contig (@Contigs){
			$TYMEFLIES_contig2MAG{$contig} = $mag;
		}
	}
}
close IN;

my %Contig2MAG = %TYMEFLIES_contig2MAG;

### Step 2.1.4 Store all contigs that are determined to be viral sequence mis-binning
my %Contig_viral_seq_in_MAGs = (); # $contig => $mag (the MAG that contains the contig)
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
	}
	close INN;
}
close IN;

### Step 2.1.5 Read hmmsearch result
my %KO2microbial_pros = (); # $ko => $microbial_pros (collection of $microbial_pro, separated by "\,")
open IN, "find Host_prediction/Find_host_based_on_AMG/ -name '*.hmmsearch_result.txt' | ";
while (<IN>){
	chomp;
	my $file = $_;
	my ($ko) = $file =~ /^.+\/(K.+?)\.hmmsearch_result\.txt/;
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
			my $microbial_pro = $tmp[0];
			my $score_full = $tmp[5];
			my $score_domain = $tmp[8];
			my ($microbial_scf) = $microbial_pro =~ /^(.+)\_\d+?$/;
			# Exclude all contigs/scaffolds that are determined to be viral sequence mis-binning
			if (! exists $Contig_viral_seq_in_MAGs{$microbial_scf}){ 
				if ($type eq "full"){
					if ($score_full >= $cutoff){
						if (!exists $KO2microbial_pros{$ko}){
							$KO2microbial_pros{$ko} = $microbial_pro;
						}else{
							$KO2microbial_pros{$ko} .= "\,".$microbial_pro;
						}
					}
				}else{
					if ($score_domain >= $cutoff){
						if (!exists $KO2microbial_pros{$ko}){
							$KO2microbial_pros{$ko} = $microbial_pro;
						}else{
							$KO2microbial_pros{$ko} .= "\,".$microbial_pro;
						}				
					}
				}
			}
		}
	}
	close INN;
}
close IN;

# Step 3 Get each KO sequences (contain both viral and microbial proteins) and write down the KO sequences
my %Microbial_protein_seq = _store_seq("/storage1/data11/TYMEFLIES_phage/Robin_MAGs/TYMEFLIES_MAGs_all_seq.faa");
my %Viral_protein_seq = _store_seq("/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/All_phage_genome.faa");

foreach my $ko (sort keys %KO2amg_pro){
	my %Seq_ko_viral_pro = (); # Store the viral protein sequences for this KO
	my %Seq_ko_microbial_pro = (); # Store the microbial protein sequences for this KO

	# Store viral protein sequences
	my @AMG_pro = split (/\,/,$KO2amg_pro{$ko});
	foreach my $amg_pro (@AMG_pro){
		my $amg_pro_new = ">$amg_pro";
		$Seq_ko_viral_pro{$amg_pro_new} = $Viral_protein_seq{$amg_pro_new};
	}
	
	# Store microbial protein sequences
	if (exists $KO2microbial_pros{$ko}){
		my @Microbial_pros = split (/\,/, $KO2microbial_pros{$ko});
		foreach my $microbial_pro (@Microbial_pros){
			my $microbial_pro_new = ">$microbial_pro";
			$Seq_ko_microbial_pro{$microbial_pro_new} = $Microbial_protein_seq{$microbial_pro_new};
		}	
	}
	
	# Write down each KO viral protein sequences
	open OUT, ">Host_prediction/Find_host_based_on_AMG/AMG_KO2all_viral_pros.$ko.faa";
	foreach my $key (sort keys %Seq_ko_viral_pro){
		print OUT "$key\n$Seq_ko_viral_pro{$key}\n";
	}
	close OUT;
	
	# Write down each KO microbial protein sequences
	if (%Seq_ko_microbial_pro){
		open OUT, ">Host_prediction/Find_host_based_on_AMG/AMG_KO2all_microbial_pros.$ko.faa";
		foreach my $key (sort keys %Seq_ko_microbial_pro){
			print OUT "$key\n$Seq_ko_microbial_pro{$key}\n";
		}
		close OUT;	
	}
}

# Step 4 Perform viral-vs-microbial blastp
open IN, "ls Host_prediction/Find_host_based_on_AMG/AMG_KO2all_microbial_pros.*.faa |";
while (<IN>){
	chomp;
	my $faa_file = $_;
	my ($faa_name) = $faa_file =~ /^.+\/(.+?)\.faa/;
	my $faa_name_viral = $faa_name; $faa_name_viral =~ s/microbial/viral/g;
	if (!(-e "Host_prediction/Find_host_based_on_AMG/$faa_name.all_vs_all_blastp_result.txt")){
		# Make diamond db
		`diamond makedb --in $faa_file --db Host_prediction/Find_host_based_on_AMG/$faa_name --threads 10`;
		# Perform all-vs-all blastp
		`diamond blastp --query Host_prediction/Find_host_based_on_AMG/$faa_name_viral.faa --db Host_prediction/Find_host_based_on_AMG/$faa_name --out Host_prediction/Find_host_based_on_AMG/$faa_name.viral_vs_microbial_blastp_result.txt --outfmt 6 --evalue 1e-5 --max-target-seqs 100 --query-cover 70 --subject-cover 70 --threads 10 --id 70 -b 16 -c 1`;
	}
}
close IN;
=cut
# Step 5 Store microbial pro to tax hash
## Step 5.1 Store the information from TYMEFLIES_all_MAGs_stat.txt, incude tax and scaffolds
my %TYMEFLIES_MAG_stat = (); # $mag => [0] GTDB tax [1] scaffolds
open IN, "/storage1/data11/TYMEFLIES_phage/Robin_MAGs/Robin_MAG_stat.CheckM_passed.txt";
while (<IN>){
	chomp;
	if (!/^IMG/){
		my @tmp = split (/\t/);
		my $mag = $tmp[0];
		my $scfs = $tmp[4];
		my $tax = $tmp[1];
		$TYMEFLIES_MAG_stat{$mag}[0] = $tax;
		$TYMEFLIES_MAG_stat{$mag}[1] = $scfs;
	}
}
close IN;

## Step 5.2 Delete MAGs that are Noncyanobacteria_but_contain_psbAD 
my %Noncyanobacteria_but_contain_psbAD = (); # $mag => 1
open IN, "/storage1/data11/TYMEFLIES_phage/Check_non-cyanobacteria/Noncyanobacteria_but_contain_psbAD.txt";
while (<IN>){
	chomp;	
	my $mag = $_;
	$Noncyanobacteria_but_contain_psbAD{$mag} = 1;
}
close IN;

## Step 5.3 Store microbial pro to tax hash
my %Microbial_pro2tax = (); # $microbial_pro => $tax
open IN, "find /storage1/data11/TYMEFLIES_phage/Robin_MAGs/All_passed_MAGs/ -name '*.faa'| grep -v 'all' |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($mag) = $file =~ /^.+\/(.+?)\.faa/;
	my $tax = $TYMEFLIES_MAG_stat{$mag}[0];
	if (! exists $Noncyanobacteria_but_contain_psbAD{$mag}){
		open INN, "$file";
		while (<INN>){
			chomp;
			if (/^>/){
				my ($microbial_pro) = $_ =~ /^>(.+?)\s/; # Only get the microbial protein before the first whitespace
				$Microbial_pro2tax{$microbial_pro} = $tax;
			}
		}
		close INN;
	}
}
close IN;

# Step 6 Parse the viral-vs-microbial blastp result
## Step 6.1 Get sequence identity_n_score hash
my %Viral_pro2microbial_pro2identity_n_score = (); # $viral_pro => $microbial_pro => $identity_n_score
my %Viral_pro2microbial_pro_all = (); # $viral_pro => $microbial_pro_all (collection of $microbial_pro, separated by "\,")
open IN, "ls Host_prediction/Find_host_based_on_AMG/*.viral_vs_microbial_blastp_result.txt | ";
while (<IN>){
	chomp;
	my $file = $_;
	if (-s "$file"){
		open INN, "$file";
		while (<INN>){
			chomp;
			my @tmp = split (/\t/);
			my $query = $tmp[0];
			my $target = $tmp[1];
			my $iden = $tmp[2];
			my $score = $tmp[-1];
			
			if ($query =~ /__vRhyme_/ and $target =~ /^Ga/){
				$Viral_pro2microbial_pro2identity_n_score{$query}{$target} = $iden."\t".$score;
				if (!exists $Viral_pro2microbial_pro_all{$query}){
					$Viral_pro2microbial_pro_all{$query} = $target;
				}else{
					$Viral_pro2microbial_pro_all{$query} .= "\,".$target;
				}
			}
		}
		close INN;
	}
}
close IN;

## Step 6.2 Get each viral pro hit (the minimal identity cutoff of viral pro to microbial pro hit is 60%)
my %Viral_pro2microbial_pro_hits = (); # $viral_pro => $microbial_pro_hits (collection of $microbial_pro_hit, separated by "\t")
foreach my $viral_pro (sort keys %Viral_pro2microbial_pro2identity_n_score){
	my @Microbial_pro_hits = (); # Only store top 10 microbial pro hits with the highest bitscores
	my %Microbial_pro_hit2score = (); # Store the bitscore for each microbial pro hit
	my @Microbial_pro_all = split (/\,/, $Viral_pro2microbial_pro_all{$viral_pro});
	foreach my $microbial_pro (@Microbial_pro_all){
		if (exists $Viral_pro2microbial_pro2identity_n_score{$viral_pro}{$microbial_pro}){
			my ($iden, $score) = $Viral_pro2microbial_pro2identity_n_score{$viral_pro}{$microbial_pro} =~ /^(.+?)\t(.+?)$/;
			if ($iden >= 60){
				$Microbial_pro_hit2score{$microbial_pro} = $score;
			}
		}
	}
	
	# Reorder the %Microbial_pro_hit2score by its values
	my @Microbial_pro_hit2score = reverse sort { $Microbial_pro_hit2score{$a} <=> $Microbial_pro_hit2score{$b} } keys %Microbial_pro_hit2score;
	if ((scalar @Microbial_pro_hit2score) <= 10){
		@Microbial_pro_hits = @Microbial_pro_hit2score;
	}else{
		splice(@Microbial_pro_hit2score, 10); # Splice the orignal array
		@Microbial_pro_hits = @Microbial_pro_hit2score;
	}
	
	if (@Microbial_pro_hits){
		$Viral_pro2microbial_pro_hits{$viral_pro} = join("\t",@Microbial_pro_hits);
	}
}

# Step 7 Get Viral_gn2host_tax hash
my %Viral_gn2viral_pro = (); # $viral_gn => $viral_pros (collection of $viral_pro, separated by "\,")
foreach my $viral_pro (sort keys %Viral_pro2microbial_pro_hits){
	my ($viral_gn) = $viral_pro =~ /^(.+?\_\_.+?)\_\_/;
	if (!exists $Viral_gn2viral_pro{$viral_gn}){
		$Viral_gn2viral_pro{$viral_gn} = $viral_pro;
	}else{
		$Viral_gn2viral_pro{$viral_gn} .= "\,".$viral_pro;
	}
}

my %Viral_gn2host_tax = (); # $viral_gn => $host_tax
foreach my $viral_gn (sort keys %Viral_gn2viral_pro){
	my $host_tax = "";
	
	my %Microbial_pro_hits = (); # $microbial_pro_hit => 1
	my @Viral_pros = split (/\,/, $Viral_gn2viral_pro{$viral_gn});
	foreach my $viral_pro (@Viral_pros){
		if (exists $Viral_pro2microbial_pro_hits{$viral_pro}){
			my @Microbial_pro_hits_this_time = split (/\t/, $Viral_pro2microbial_pro_hits{$viral_pro});
			if (@Microbial_pro_hits_this_time){
				foreach my $key (@Microbial_pro_hits_this_time){
					$Microbial_pro_hits{$key} = 1;
				}
			}
		}
	}
	
	my @Taxs = (); # Store all the tax 
	if (%Microbial_pro_hits){
		foreach my $microbial_pro_hit (sort keys %Microbial_pro_hits){
			if (exists $Microbial_pro2tax{$microbial_pro_hit}){
				my $tax = $Microbial_pro2tax{$microbial_pro_hit};			
				push @Taxs, $tax;
			}
		}
	}

	if (@Taxs){
		$host_tax = _find_consensus_lineages_based_on_each_rank(@Taxs);
	}

	if ($host_tax){
		$Viral_gn2host_tax{$viral_gn} = $host_tax;
	}
}

# Step 8 Write down Viral_gn2host_tax hash
open OUT, ">Host_prediction/Phage_gn2host_tax_based_on_AMG.txt";
foreach my $gn (sort keys %Viral_gn2host_tax){
	print OUT "$gn\t$Viral_gn2host_tax{$gn}\n";
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


