#!/usr/bin/perl

use strict;
use warnings;

# AIM: Grep all AMG proteins and run against dbCAN2 and MEROPS databases to get annotated

# Step 1. Grep all AMG proteins
my %AMG_pro2ko = (); # $amg_pro => $ko
my %AMG_pro2line = (); # $amg_pro => $line
my $header = ""; # Store the header of AMG_summary.txt
open IN, "AMG_analysis/AMG_summary.txt";
while (<IN>){
	chomp;
	if (/^Pro/){
		$header = $_;
	}else{	
		my $line = $_;
		my @tmp = split (/\t/, $line);
		my $amg_pro = $tmp[0];
		my $ko = $tmp[2];
		$AMG_pro2ko{$amg_pro} = $ko;
		$AMG_pro2line{$amg_pro} = $line;
	}
}
close IN;
=pod
# Store all phage faa and write down All_AMG_pro.faa
`find /storage1/data11/TYMEFLIES_phage/*/vRhyme_best_bins_fasta_parsed -name '*.faa' -exec cat {} + > /storage1/data11/TYMEFLIES_phage/All_phage_pro.faa`;
my %All_phage_pro = _store_seq("/storage1/data11/TYMEFLIES_phage/All_phage_pro.faa");
`rm /storage1/data11/TYMEFLIES_phage/All_phage_pro.faa`;

open OUT, ">All_AMG_pro.faa";
foreach my $key (sort keys %All_phage_pro){
	my ($key_clean) = $key =~ /^>(.+?)$/;
	if (exists $AMG_pro2ko{$key_clean}){
		print OUT "$key\n$All_phage_pro{$key}\n";
	}
}
close OUT;

# Step 2. Run dbCAN2 to find CAZymes
`mkdir AMG_analysis/dbCAN2_and_MEROPS_analysis`;

# Split All_AMG_pro.faa into mutiple faa files to allow parallel run
`mkdir AMG_analysis/dbCAN2_and_MEROPS_analysis/split_fsa`;
`perl /storage1/data11/TYMEFLIES_phage/split_multifasta.pl --in All_AMG_pro.faa --output_dir=AMG_analysis/dbCAN2_and_MEROPS_analysis/split_fsa --seqs_per_file=1000`;
`rm All_AMG_pro.faa`;

open OUT, ">tmp.run_dbCAN2_in_parallel.sh";
open IN, "ls AMG_analysis/dbCAN2_and_MEROPS_analysis/split_fsa/*.fsa |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($file_name) = $file =~ /split_fsa\/(.+?)\.fsa/;
	print OUT "hmmscan --domtblout AMG_analysis/dbCAN2_and_MEROPS_analysis/$file_name.dbCAN2.out.dm --cpu 1 /slowdata/data1/Genome_profile_software/dbCAN2/dbCAN-fam-HMMs.txt $file > AMG_analysis/dbCAN2_and_MEROPS_analysis/$file_name.dbCAN2.out;";
	print OUT "python /slowdata/data1/Genome_profile_software/Accessory_scripts/hmmscan-parser-dbCANmeta.py AMG_analysis/dbCAN2_and_MEROPS_analysis/$file_name.dbCAN2.out.dm > AMG_analysis/dbCAN2_and_MEROPS_analysis/$file_name.dbCAN2.out.dm.ps\n";
}
close IN;
close OUT;

`cat tmp.run_dbCAN2_in_parallel.sh | parallel -j 5`;
`rm tmp.run_dbCAN2_in_parallel.sh`;
=cut
# Step 3. Parse dbCAN result
## Step 3.1 Store CAZy_map
my %CAZy_map = (); # $cazyme_id => [0] $cazyme_enzyme [1] $cazyme_substrate [2] $cazyme_class_of_substrate
open IN, "/slowdata/data1/Genome_profile_software/METABOLIC_template_and_database/CAZy_map.txt";
while (<IN>){
	chomp;
	if (!/^Family/){
		my @tmp = split (/\t/);
		my $cazyme_id = $tmp[0];
		my $cazyme_enzyme = $tmp[1];
		my $cazyme_substrate = $tmp[2];
		my $cazyme_class_of_substrate = $tmp[3];
		$CAZy_map{$cazyme_id}[0] = $cazyme_enzyme;
		$CAZy_map{$cazyme_id}[1] = $cazyme_substrate;
		$CAZy_map{$cazyme_id}[2] = $cazyme_class_of_substrate;
	}
}
close IN;

my %AMG_pro2cazyme = (); # $amg_pro => [0] $cazyme_id; if there multiple $cazyme_id hits for one protein, then store all, and separate them by "|"
                                      #[1] $cazyme_detail; if there multiple $cazyme_detail hits for one protein, then store all, and separate them by "|"
                                      # If the corresponding $cazyme_detail could be found then write "NA"
open IN, "cat AMG_analysis/dbCAN2_and_MEROPS_analysis/*.dbCAN2.out.dm.ps |";
while (<IN>){
	chomp;
	if (/^GH|^PL/){
		my @tmp = split (/\t/);
		my $amg_pro = $tmp[2];
		my ($cazyme_id) = $tmp[0] =~ /^(.+?)\.hmm/;
		# Change GH13_28 to GH130_28
		my $cazyme_id_new = "";
		my $cazyme_id_front = ""; # Store only the cazyme id if in front of "_"
		if ($cazyme_id =~ /_/){
			my ($cazyme_id_p1,$cazyme_id_p2,$cazyme_id_p3) = $cazyme_id =~ /^(\D+?)(\d+?)\_(\d+?)$/;
			$cazyme_id_p2 = (sprintf "%03d", $cazyme_id_p2); # fill numbers into 3-position containing numbers by adding "0" in the front
			$cazyme_id_new = $cazyme_id_p1.$cazyme_id_p2."_".$cazyme_id_p3;
			$cazyme_id_front = $cazyme_id_p1.$cazyme_id_p2;
		}else{
			my ($cazyme_id_p1,$cazyme_id_p2) = $cazyme_id =~ /^(\D+?)(\d+?)$/;
			$cazyme_id_p2 = (sprintf "%03d", $cazyme_id_p2); # fill numbers into 3-position containing numbers by adding "0" in the front
			$cazyme_id_new = $cazyme_id_p1.$cazyme_id_p2;
			$cazyme_id_front = $cazyme_id_p1.$cazyme_id_p2;
		}
		
		if (!exists $AMG_pro2cazyme{$amg_pro}[0]){
			$AMG_pro2cazyme{$amg_pro}[0] = $cazyme_id_new;
		}else{
			$AMG_pro2cazyme{$amg_pro}[0] .= "\|".$cazyme_id_new;
		}
		
		my $cazyme_detail = "NA\tNA\tNA";
		if (exists $CAZy_map{$cazyme_id_front}[0]){
			$cazyme_detail = $CAZy_map{$cazyme_id_front}[0]."\t".$CAZy_map{$cazyme_id_front}[1]."\t".$CAZy_map{$cazyme_id_front}[2];
		}
		
		if (!exists $AMG_pro2cazyme{$amg_pro}[1]){
			$AMG_pro2cazyme{$amg_pro}[1] = $cazyme_detail;
		}else{
			$AMG_pro2cazyme{$amg_pro}[1] .= "\|".$cazyme_detail;
		}		
	}
}
close IN;
=pod
# Step 4. Run MEROPS to find peptidases
open OUT, ">tmp.run_MEROPS_in_parallel.sh";
open IN, "ls AMG_analysis/dbCAN2_and_MEROPS_analysis/split_fsa/*.fsa |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($file_name) = $file =~ /split_fsa\/(.+?)\.fsa/;
	print OUT "diamond blastp -d /slowdata/data1/Genome_profile_software/MEROPS/pepunit.db -q $file -o AMG_analysis/dbCAN2_and_MEROPS_analysis/$file_name.MEROPSout.m8 -k 1 -e 1e-10 --query-cover 80 --id 50 --quiet -p 1 > /dev/null\n";
}
close IN;
close OUT;

`cat tmp.run_MEROPS_in_parallel.sh | parallel -j 5`;
`rm tmp.run_MEROPS_in_parallel.sh`;
=cut
# Step 5. Parse MEROPS result
## Step 5.1 Store MEROPS map 
my %MEROPS_map; # MER id => all line; for example: MER0000001 => all line
open IN, "/slowdata/data1/Genome_profile_software/MEROPS/pepunit.lib";
while (<IN>){
	chomp;
	if (/>/){
		$_ =~ tr/\015//d;
		my ($mer_id) = $_ =~ /^>(.+?)\s/;
		$MEROPS_map{$mer_id} = $_; 
	}
}
close IN;

my %AMG_pro2merops = (); # $amg_pro => [0]$merops_id; if there multiple $merops_id hits for one protein, then store all of them, and separate them by "|"
                                      #[1]$merops_detail; if there multiple $merops_detail hits for one protein, then store all of them, and separate them by "|"
open IN, "cat AMG_analysis/dbCAN2_and_MEROPS_analysis/*.m8 |";
while (<IN>){
	chomp;
	my @tmp = split(/\t/);
	my ($merops_detail,$merops_id) = $MEROPS_map{$tmp[1]} =~ /\s\-\s(.+?)\#(.+?)\#/; 
	my $amg_pro = $tmp[0];
	if (!exists $AMG_pro2merops{$amg_pro}[0]){
		$AMG_pro2merops{$amg_pro}[0] = $merops_id;
	}else{
		$AMG_pro2merops{$amg_pro}[0] .= "\|".$merops_id;
	}

	if (!exists $AMG_pro2merops{$amg_pro}[1]){
		$AMG_pro2merops{$amg_pro}[1] = $merops_detail;
	}else{
		$AMG_pro2merops{$amg_pro}[1] .= "\|".$merops_detail;
	}	
}
close IN;

# Step 6. Add dbCAN2 and MEROPS annotation results to the AMG summary table
my %AMG_pro2info = (); # $amg_pro => $line (new, add dbCAN2 and MEROPS results)
                       # Add 1) $cazyme_id 2) $cazyme_detail (contains 3 columns) 3) $merops_id 4) $merops_detail

foreach my $amg_pro (sort keys %AMG_pro2line){
	my $cazyme_id = $AMG_pro2cazyme{$amg_pro}[0];
	my $cazyme_detail = $AMG_pro2cazyme{$amg_pro}[1];
	my $merops_id = $AMG_pro2merops{$amg_pro}[0];
	my $merops_detail = $AMG_pro2merops{$amg_pro}[1];
	$AMG_pro2info{$amg_pro}[0] = $cazyme_id;
	$AMG_pro2info{$amg_pro}[1] = $cazyme_detail;
	$AMG_pro2info{$amg_pro}[2] = $merops_id;
	$AMG_pro2info{$amg_pro}[3] = $merops_detail;
}

# Step 7. Write down new AMG_summary.txt
$header = $header."\t"."CAZyme ID"."\t"."Enzyme (CAZyme)"."\t"."Substrate (CAZyme)"."\t"."Class of Substrate (CAZyme)"."\t"."MEROPS ID"."\t"."MEROPS details";
open OUT, ">AMG_analysis/AMG_summary_new.txt";
print OUT "$header\n";
foreach my $amg_pro (sort keys %AMG_pro2info){
	my $line = $AMG_pro2line{$amg_pro};
	my $cazyme_id = ($AMG_pro2info{$amg_pro}[0])? $AMG_pro2info{$amg_pro}[0]:"NA";
	my $cazyme_detail = ($AMG_pro2info{$amg_pro}[1])? $AMG_pro2info{$amg_pro}[1]:"NA\tNA\tNA";
	my $merops_id = ($AMG_pro2info{$amg_pro}[2])? $AMG_pro2info{$amg_pro}[2]:"NA";
	my $merops_detail = ($AMG_pro2info{$amg_pro}[3])? $AMG_pro2info{$amg_pro}[3]:"NA";	
	
	print OUT "$line\t$cazyme_id\t$cazyme_detail\t$merops_id\t$merops_detail\n";
}
close OUT;

# Step 8. Replace the old one
## Sort the result by date
`cat AMG_analysis/AMG_summary_new.txt | sort -k 2 -n > tmp`;
`mv tmp AMG_analysis/AMG_summary_new.txt`;

## Replace the old summary file
`rm AMG_analysis/AMG_summary.txt`;
`mv AMG_analysis/AMG_summary_new.txt AMG_analysis/AMG_summary.txt`;


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
