#!/usr/bin/perl

use strict;
use warnings;

# AIM: Get Viral_gn_with_psbAD2host_mag.by_AMG_match result
#      viral_gn_with_psbAD => host_mags by the method of AMG match


# Step 1 Store microbial pro to mag hash
## Step 1.1 Store the information from TYMEFLIES_all_MAGs_stat.txt, incude tax and scaffolds
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

## Step 1.2 Store microbial pro to mag hash
my %Microbial_pro2mag = (); # $microbial_pro => $mag
open IN, "find /storage1/data11/TYMEFLIES_phage/Robin_MAGs/All_passed_MAGs/ -name '*.faa'| grep -v 'all' |";
while (<IN>){
	chomp;
	my $file = $_;
	my ($mag) = $file =~ /^.+\/(.+?)\.faa/;
	open INN, "$file";
	while (<INN>){
		chomp;
		if (/^>/){
			my ($microbial_pro) = $_ =~ /^>(.+?)\s/; # Only get the microbial protein before the first whitespace
			$Microbial_pro2mag{$microbial_pro} = $mag;
		}
	}
	close INN;
}
close IN;

# Step 2 Parse the viral-vs-microbial blastp result
## Step 2.1 Get sequence identity_n_score hash
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

## Step 2.2 Get each viral pro hit (the minimal identity cutoff of viral pro to microbial pro hit is 90%)
my %Viral_pro2microbial_pro_hits = (); # $viral_pro => $microbial_pro_hits (collection of $microbial_pro_hit, separated by "\t")
foreach my $viral_pro (sort keys %Viral_pro2microbial_pro2identity_n_score){
	my @Microbial_pro_hits = (); # Only store top 20 microbial pro hits with the highest bitscores
	my %Microbial_pro_hit2score = (); # Store the bitscore for each microbial pro hit
	my @Microbial_pro_all = split (/\,/, $Viral_pro2microbial_pro_all{$viral_pro});
	foreach my $microbial_pro (@Microbial_pro_all){
		if (exists $Viral_pro2microbial_pro2identity_n_score{$viral_pro}{$microbial_pro}){
			my ($iden, $score) = $Viral_pro2microbial_pro2identity_n_score{$viral_pro}{$microbial_pro} =~ /^(.+?)\t(.+?)$/;
			if ($iden >= 90){
				$Microbial_pro_hit2score{$microbial_pro} = $score;
			}
		}
	}
	
	# Reorder the %Microbial_pro_hit2score by its values
	my @Microbial_pro_hit2score = reverse sort { $Microbial_pro_hit2score{$a} <=> $Microbial_pro_hit2score{$b} } keys %Microbial_pro_hit2score;
	if ((scalar @Microbial_pro_hit2score) <= 20){
		@Microbial_pro_hits = @Microbial_pro_hit2score;
	}else{
		splice(@Microbial_pro_hit2score, 20); # Splice the orignal array
		@Microbial_pro_hits = @Microbial_pro_hit2score;
	}
	
	if (@Microbial_pro_hits){
		$Viral_pro2microbial_pro_hits{$viral_pro} = join("\t",@Microbial_pro_hits);
	}
}

# Step 3 Get Viral_gn2host_mag hash
my %Viral_gn2viral_pro = (); # $viral_gn => $viral_pros (collection of $viral_pro, separated by "\,")
foreach my $viral_pro (sort keys %Viral_pro2microbial_pro_hits){
	my ($viral_gn) = $viral_pro =~ /^(.+?\_\_.+?)\_\_/;
	if (!exists $Viral_gn2viral_pro{$viral_gn}){
		$Viral_gn2viral_pro{$viral_gn} = $viral_pro;
	}else{
		$Viral_gn2viral_pro{$viral_gn} .= "\,".$viral_pro;
	}
}

my %Viral_gn2host_mag = (); # $viral_gn => $mags (collection of $mag, separated by "\,")
foreach my $viral_gn (sort keys %Viral_gn2viral_pro){	
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
	
	my @MAGs = (); # Store all the mags 
	if (%Microbial_pro_hits){
		foreach my $microbial_pro_hit (sort keys %Microbial_pro_hits){
			my $mag = $Microbial_pro2mag{$microbial_pro_hit};			
			push @MAGs, $mag;
		}
	}

	if (@MAGs){
		$Viral_gn2host_mag{$viral_gn} = join("\,", @MAGs);
	}
}

# Step 4 Find non-cyanobacteria MAG hosts that are corresponding to 
#        viral genome with psbAD 

## Store Viral_gn_with_psbAD_connected_to_noncyanobacteria_host.txt
my %Viral_gn_with_psbAD_connected_to_noncyanobacteria_host = ();
open IN, "Viral_gn_with_psbAD_connected_to_noncyanobacteria_host.txt";
while (<IN>){
	chomp;
	my $viral_gn = $_;
	$Viral_gn_with_psbAD_connected_to_noncyanobacteria_host{$viral_gn} = 1;
}
close IN;

my %Viral_gn_with_psbAD2host_mag = (); # $viral_gn (with psbAD inside) => $mags (collection of $mag separated by "\,", the host) and $taxs (collection of $tax separated by "\,") 
foreach my $viral_gn (sort keys %Viral_gn_with_psbAD_connected_to_noncyanobacteria_host){
	my %MAGs_all = (); # Store all the MAGs 
	if (exists $Viral_gn2host_mag{$viral_gn}){
		my $mags = $Viral_gn2host_mag{$viral_gn};
		my @MAGs = split(/\,/, $mags);
		foreach my $mag (@MAGs){
			$MAGs_all{$mag} = 1;
		}
		
	}
	
	if (%MAGs_all){
		my @MAGs = sort keys %MAGs_all;
		my $mags = join("\,", @MAGs);
		my @Taxs = ();
		foreach my $mag (@MAGs){
			my $tax = $TYMEFLIES_MAG_stat{$mag}[0];
			push @Taxs, $tax;
		}
		my $taxs = join("\,", @Taxs);
		
		$Viral_gn_with_psbAD2host_mag{$viral_gn} = $mags."\|".$taxs;
	}	
}

## Write down the result
open OUT, ">Viral_gn_with_psbAD2host_mag.by_AMG_match.txt";
foreach my $viral_gn (sort keys %Viral_gn_with_psbAD2host_mag){
	print OUT "$viral_gn\t$Viral_gn_with_psbAD2host_mag{$viral_gn}\n";
}
close OUT;
