#!/usr/bin/perl

use strict;
use warnings;

############################################################################################
# To get new vRhyme bins, due to that                                                      #

# 1. Prophage identified by VIBRANT should be excluded from binning                        #

# 2. Two or more lysogenic (non-prophage) phage scaffolds can not be in the same bin       #                      

# 3. Phage scaffolds identified by CheckV as "Complete" should be excluded from binning    #

# 4. The maximum number of bin redundancy should be <= 1                                   #
############################################################################################


## Step 1: Make new membership tsv
## Step 2: Make new bins

## Step 1: Make new membership tsv
my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

foreach my $img_id (sort keys %IMGID){
	
	# Get all scaffold state
	my %Scf2state = (); # $scf => $state (either "prophage", "lytic", or "lysogenic_but_non-prophage"); $scf here contains fragment inside
	open IN, "$img_id/VIBRANT_$img_id.a.v2.min2000/phage_results.min2000.txt";
	while (<IN>){
		chomp;
		my @tmp = split (/\t/);
		my $scf = $tmp[0];
		my $state = $tmp[1];
		$Scf2state{$scf} = $state;
	}
	close IN;
	
	# Store membership
	my %Membership = ();  # $scf => [0] $bin [1] $state [2] $note [3] $redundancy
	my $membership_file_name = ""; # For example, "vRhyme_best_bins.0.membership"
	my %Bin = (); # $bin => $scf separated by "\t"; $scf here contains fragment inside
	open IN, "ls $img_id/vRhyme_result/*membership.tsv |";
	while (<IN>){
		chomp;
		my $file = $_;
		($membership_file_name) = $file =~ /vRhyme\_result\/(.+?)\.tsv/;
		open INN, "$file";
		while (<INN>){
			chomp;
			if (!/^scaffold/){
				my @tmp = split (/\t/);
				my $scf = $tmp[0];
				my $bin = $tmp[1];
				if (!exists $Bin{$bin}){
					$Bin{$bin} = $scf;
				}else{
					$Bin{$bin} .= "\t".$scf;
				}
				$Membership{$scf}[0] = $bin;
				$Membership{$scf}[1] = $Scf2state{$scf};
			}
		}
		close INN;
	}
	close IN;
	
	# Store CheckV "Complete" info
	my %Scf2checkV_complete = (); # $scf => 1; $scf here contains fragment inside
	open IN, "$img_id/CheckV_phage_scaffold/quality_summary.tsv";
	while (<IN>){
		chomp;
		if (!/^contig_id/){
			my @tmp = split (/\t/);
			if ($tmp[7] eq "Complete"){
				$Scf2checkV_complete{$tmp[7]} = 1;
			}
		}
	}
	close IN;
	
	# Store redundancy
	my %Bin2redundancy = ();  # $bin => $redundancy
	open IN, "ls $img_id/vRhyme_result/*summary.tsv |";
	while (<IN>){
		chomp;
		my $file = $_;
		open INN, "$file";
		while (<INN>){
			chomp;
			if (!/^bin/){
				my @tmp = split (/\t/);
				my $bin = $tmp[0];
				my $redundancy = $tmp[3];
				$Bin2redundancy{$bin} = $redundancy;
			}
		}
		close INN;
	}
	close IN;
	
	foreach my $bin (sort keys %Bin){
		my @Scfs = split (/\t/, $Bin{$bin});
		
		# Split bins with >= 1 prophage(s)
		my $prophage_num = 0;	
		foreach my $scf (@Scfs){
			if ($Membership{$scf}[1] eq "prophage"){
				$prophage_num++;
			}
		}		
		
		if ($prophage_num >= 1){
			foreach my $scf (@Scfs){
				$Membership{$scf}[2] = "split; contains >=1 prophage(s)";
			}
		}		
		
		# Split bins with >= 2 lysogenic (non-prophage) phage scaffolds
		my $lysogenic_num = 0;
		foreach my $scf (@Scfs){
			if ($Membership{$scf}[1] eq "lysogenic_but_non-prophage"){
				$lysogenic_num++;
			}
		}
		
		if ($lysogenic_num >= 2){
			foreach my $scf (@Scfs){
				if ($Membership{$scf}[2]){
					$Membership{$scf}[2] .= " | split; contains >=2 lysogenic but non-prophage scaffolds";
				}else{
					$Membership{$scf}[2] = "split; contains >=2 lysogenic but non-prophage scaffolds";
				}
			}
		}

		# Split bins containing phage scaffolds identified by CheckV as "Complete"
		my $checkV_complete_phage_scaffold_num = 0;
		foreach my $scf (@Scfs){
			if (exists $Scf2checkV_complete{$scf}){
				$checkV_complete_phage_scaffold_num++;
			}
		}

		if ($checkV_complete_phage_scaffold_num >= 1){
			foreach my $scf (@Scfs){
				if ($Membership{$scf}[2]){
					$Membership{$scf}[2] .= " | split; contains >= 1 CheckV Complete phage(s)";
				}else{
					$Membership{$scf}[2] = "split; contains >= 1 CheckV Complete phage(s)";
				}
			}
		}			

		# The maximum redundancy should be <= 1
		my $bin_redundancy = $Bin2redundancy{$bin};
		if ($bin_redundancy > 1){
			foreach my $scf (@Scfs){
				if ($Membership{$scf}[2]){
					$Membership{$scf}[2] .= " | split; bin redundancy > 1";
				}else{
					$Membership{$scf}[2] = "split; bin redundancy > 1";
				}
			}
		}
	}
	
	# Write new membership file
	open OUT, ">$img_id/vRhyme_result/$membership_file_name.parsed.tsv";
	print OUT "IMG ID\tscaffold\tbin\tstate\tnote\tredundancy\n";
	foreach my $scf (sort keys %Membership){
		my $bin = $Membership{$scf}[0];
		my $state = $Membership{$scf}[1];
		my $note = "n/a";
		if ($Membership{$scf}[2]){
			$note = $Membership{$scf}[2];
		}
		my $redundancy = $Bin2redundancy{$bin};
		
		print OUT "$img_id\t$scf\t$bin\t$state\t$note\t$redundancy\n";
	}
	close OUT;
	
	# Sort the resulted membership file (a.k.a, $membership_file_name.parsed.tsv)
	`cat $img_id/vRhyme_result/$membership_file_name.parsed.tsv | sort  -k3,3n > $img_id/vRhyme_result/tmp`;
	`mv $img_id/vRhyme_result/tmp  $img_id/vRhyme_result/$membership_file_name.parsed.tsv`;
	
	
	# Step 2: Make new bins and unbinned scaffold
	# Here new bin format was made like this:
	# For fasta file: sequence head is like ">3300020489__vRhyme_1__Ga0207910_100402"
	# For faa file: sequence head is like ">3300020489__vRhyme_1__Ga0207910_100402_1"
	# For ffn file: sequence head is like ">3300020489__vRhyme_1__Ga0207910_100402_1"
	# Unbinned scaffold format was made like this:
	# For fasta file: sequence head is like ">3300020489__vRhyme_unbinned__Ga0207910_100402"
	# For faa file: sequence head is like ">3300020489__vRhyme_unbinned__Ga0207910_100402_1"
	# For ffn file: sequence head is like ">3300020489__vRhyme_unbinned__Ga0207910_100402_1"
	my %Seq_phage_fna = _store_seq("$img_id/VIBRANT_$img_id.a.v2.min2000/$img_id.a.phage_scaffold.fna");
	my %Seq_phage_faa = _store_seq("$img_id/VIBRANT_$img_id.a.v2.min2000/$img_id.a.phage_scaffold.faa");
	my %Seq_phage_ffn = _store_seq("$img_id/VIBRANT_$img_id.a.v2.min2000/$img_id.a.phage_scaffold.ffn");
	
	my %Each_bin_info = (); # Record the info of bin and unbinned scaffold (treated as a bin)
	# $each_bin => [0] $bin_scaffold_number [1] $bin_scaffold_length [2] $bin_protein_number
	# Here $each_bin is like 3300020489__vRhyme_1 and 3300020489__vRhyme_unbinned1
	
	`mkdir $img_id/vRhyme_best_bins_fasta_parsed`;
	my $bin_folder = "$img_id/vRhyme_best_bins_fasta_parsed";
	
	foreach my $bin (sort keys %Bin){ # Right bracket for each bin
		my @Scfs = split (/\t/, $Bin{$bin}); # $scf here contains "fragment" inside
		
		my $action = "keep"; # The action of to split or to keep
		foreach my $scf (@Scfs){ # $scf here contains "fragment" inside
			my $note = "n/a";
			if ($Membership{$scf}[2]){
				$note = $Membership{$scf}[2];
			}
			if ($note =~ /split/){
				$action = "split";
			}
		}
		
		if ($action eq "keep"){ # Keep the bin, change the sequence head, and write down the bin
			my %Each_bin_seq_fna = (); # Temporary hash to store each bin's fna sequences
			my %Each_bin_seq_faa = (); # Temporary hash to store each bin's faa sequences
			my %Each_bin_seq_ffn = (); # Temporary hash to store each bin's ffn sequences
			
			# Make each bin fna hash
			foreach my $scf (@Scfs){ 
				my $head_old = ">$scf"; # Old head in %Seq_phage_fna
				my $head_new = ">$img_id\_\_vRhyme\_$bin\_\_$scf"; # New head in %Each_bin_seq_fna
				$Each_bin_seq_fna{$head_new} = $Seq_phage_fna{$head_old};
				# Delete binned scaffold
				delete $Seq_phage_fna{$head_old};
			}
			my @Each_bin_seq_fna_stat = _get_stat_from_seq_hash(%Each_bin_seq_fna);
			$Each_bin_info{"$img_id\_\_vRhyme\_$bin"}[0] = $Each_bin_seq_fna_stat[0]; # Store scaffold num
			$Each_bin_info{"$img_id\_\_vRhyme\_$bin"}[1] = $Each_bin_seq_fna_stat[1]; # Store scaffold length
			
			# Make each bin faa hash
			my %Pro_head_binned = (); # The sequence head of proteins in phage bins
			foreach my $scf (@Scfs){ 
			    foreach my $head (sort keys %Seq_phage_faa){
					if ($head =~ />$scf\_/){
						$Pro_head_binned{$head} = 1;
					}
				}
			}

			foreach my $head (sort keys %Pro_head_binned){
				my $head_old = $head; my $head_old_clean = $head_old; $head_old_clean =~ s/^>//g;
				my $head_new = ">$img_id\_\_vRhyme\_$bin\_\_$head_old_clean"; # New head in %Each_bin_seq_faa
				$Each_bin_seq_faa{$head_new} = $Seq_phage_faa{$head_old};
				# Delete proteins in binned scaffold
				delete $Seq_phage_faa{$head_old};
			}
			
			my @Each_bin_seq_faa_stat = _get_stat_from_seq_hash(%Each_bin_seq_faa);
			$Each_bin_info{"$img_id\_\_vRhyme\_$bin"}[2] = $Each_bin_seq_faa_stat[0]; # Store protein num
			
			# Make each bin ffn hash
			my %Gene_head_binned = (); # The sequence head of proteins in phage bins
			%Gene_head_binned = %Pro_head_binned; # Gene head are the same as Protein head

			foreach my $head (sort keys %Gene_head_binned){
				my $head_old = $head; my $head_old_clean = $head_old; $head_old_clean =~ s/^>//g;
				my $head_new = ">$img_id\_\_vRhyme\_$bin\_\_$head_old_clean"; # New head in %Each_bin_seq_ffn
				$Each_bin_seq_ffn{$head_new} = $Seq_phage_ffn{$head_old};
				# Delete genes in binned scaffold
				delete $Seq_phage_ffn{$head_old};
			}

			# Write down fna, faa, and ffn of each bin
			open OUT, ">$bin_folder/$img_id\_\_vRhyme_$bin.fasta";
			foreach my $key (sort keys %Each_bin_seq_fna){
				print OUT "$key\n$Each_bin_seq_fna{$key}\n";
			}
			close OUT;
			
			open OUT, ">$bin_folder/$img_id\_\_vRhyme_$bin.faa";
			foreach my $key (sort keys %Each_bin_seq_faa){
				print OUT "$key\n$Each_bin_seq_faa{$key}\n";
			}			
			close OUT;
			
			open OUT, ">$bin_folder/$img_id\_\_vRhyme_$bin.ffn";
			foreach my $key (sort keys %Each_bin_seq_ffn){
				print OUT "$key\n$Each_bin_seq_ffn{$key}\n";
			}			
			close OUT;				
		}			
	} # Left bracket for each bin

	# Write down fna, faa, and ffn for unbinned scaffolds
	my $i = 1; # The unbinned order started from 1
	my %Unbinned2scf = (); # The map of unbinned order (full) to scf
	# For example: $img_id\_\_vRhyme\_unbinned1 => $scf (can contain "fragment" here)
	
	foreach my $key (sort keys %Seq_phage_fna){
		my $key_new = $key; my $key_clean = $key; $key_clean =~ s/^>//g;
		$key_new =~ s/^>/>$img_id\_\_vRhyme\_unbinned$i\_\_/g; # Change the head
		
		my $unbinned = "$img_id\_\_vRhyme\_unbinned$i"; # Store the unbinned order (full)
		$Unbinned2scf{$unbinned} = $key_clean; # Store the map of unbinned order (full) to scf
		
		my %Each_unbinned_seq_fna = (); # Store the hash for unbinned seq fna
		$Each_unbinned_seq_fna{$key_new} = $Seq_phage_fna{$key};
		
		# Store the each unbinned info
		my @Each_unbinned_seq_fna_stat = _get_stat_from_seq_hash(%Each_unbinned_seq_fna);
		$Each_bin_info{"$img_id\_\_vRhyme\_unbinned$i"}[0] = $Each_unbinned_seq_fna_stat[0]; # Store scaffold num
		$Each_bin_info{"$img_id\_\_vRhyme\_unbinned$i"}[1] = $Each_unbinned_seq_fna_stat[1]; # Store scaffold length
		
		# Write down unbinned seq fna
		open OUT, ">$bin_folder/$img_id\_\_vRhyme_unbinned$i.fasta";
		foreach my $head (sort keys %Each_unbinned_seq_fna){
			print OUT "$head\n$Each_unbinned_seq_fna{$head}\n";
		}
		close OUT;
		
		$i++; 
	}
	
	foreach my $unbinned (sort keys %Unbinned2scf){ # Analyze each unbinned scaffold
		# Store the hash for unbinned seq faa
		my %Each_unbinned_seq_faa = (); 
		my $scf = $Unbinned2scf{$unbinned};
		
		foreach my $key (sort keys %Seq_phage_faa){
			my $key_clean = $key; $key_clean =~ s/^>//g;
			if ($key_clean =~ /$scf\_/){
				my $key_new = ">$unbinned\_\_$key_clean";
				$Each_unbinned_seq_faa{$key_new} = $Seq_phage_faa{$key};
			}
		}
		
		# Store the each unbinned info
		my @Each_unbinned_seq_faa_stat = _get_stat_from_seq_hash(%Each_unbinned_seq_faa);
		$Each_bin_info{$unbinned}[2] = $Each_unbinned_seq_faa_stat[0]; # Store protein num
			
		# Write down unbinned seq faa
		open OUT, ">$bin_folder/$unbinned.faa";
		foreach my $head (sort keys %Each_unbinned_seq_faa){
			print OUT "$head\n$Each_unbinned_seq_faa{$head}\n";
		}
		close OUT;		
		
	# Store the hash for unbinned seq ffn
		my %Each_unbinned_seq_ffn = (); 
		
		foreach my $key (sort keys %Seq_phage_ffn){
			my $key_clean = $key; $key_clean =~ s/^>//g;
			if ($key_clean =~ /$scf\_/){
				my $key_new = ">$unbinned\_\_$key_clean";
				$Each_unbinned_seq_ffn{$key_new} = $Seq_phage_ffn{$key};
			}
		}

		# Write down unbinned seq ffn
		open OUT, ">$bin_folder/$unbinned.ffn";
		foreach my $head (sort keys %Each_unbinned_seq_ffn){
			print OUT "$head\n$Each_unbinned_seq_ffn{$head}\n";
		}
		close OUT;					
	}

	# Write down %Each_bin_info hash
	open OUT, ">$bin_folder/Each_bin_info.txt";
	print OUT "bin\tbin scaffold number\tbin scaffold length\tbin protein number\n";
	foreach my $key (sort keys %Each_bin_info){
		print OUT "$key\t$Each_bin_info{$key}[0]\t$Each_bin_info{$key}[1]\t$Each_bin_info{$key}[2]\n";
	}
	close OUT;
}



# subroutine

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

sub _get_stat_from_seq_hash{
	my (%Seq) = @_;
	my $seq_num = 0;
	my $seq_size = 0;
	foreach my $key (sort keys %Seq){
		$seq_num++;
		my $len = length($Seq{$key}); # The length of the sequence
		$seq_size += $len;
	}
	
	my @Return = ();
	push @Return, $seq_num;
	push @Return, $seq_size;
	
	return @Return;
}