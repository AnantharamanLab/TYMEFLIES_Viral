#!/usr/bin/perl

use strict;
use warnings;

#############################################################################################################
# Cutoff for two prophages in the same saffold to be combined and this scaffold is a entire phage scaffold: #
# 1. protein distance <= 20 proteins or nt distance <= 20 kb                                                #
# 2. a fragment gap in between both prophages (for example, prophages are fragment_1 and fragment_3)        #
# 3. host region (front, end, and gap regions counted together) <= 30%                                      #

# Cutoff for three or more prophages:                                                                       #
# 1. protein distance <= 20 proteins or nt distance <= 20 kb                                                #
# 2. host region (front, end, and gap regions counted together) <= 30%                                      #
#############################################################################################################

my $min_length = 2000; # The min length of phage scaffolds

my %IMGID = (); # $img_id => 1
open IN, "ls -l | sed -r s'/ +/\t/g' | grep ^d | cut -f 9 | grep '33' |";
while (<IN>){
	chomp;
	my $img_id = $_;
	$IMGID{$img_id} = 1;
}
close IN;

# TEST script
# Run 1 metagenome to have a test 
# %IMGID = (); $IMGID{"3300042291"} = 1;

foreach my $img_id (sort keys %IMGID){
	# Get all scaffold state
	my %Scf2state = (); # $scf => $state (either "prophage", "lytic", or "lysogenic_but_non-prophage")
	my %Seq_head_lytic = _store_seq_head_to_state("$img_id/VIBRANT_$img_id.a/VIBRANT_phages_$img_id.a/$img_id.a.phages_lytic.fna","lytic");
	my %Seq_head_lysogenic = _store_seq_head_to_state("$img_id/VIBRANT_$img_id.a/VIBRANT_phages_$img_id.a/$img_id.a.phages_lysogenic.fna","lysogenic_but_non-prophage");	
	my %Seq_head_prophage = ();
	foreach my $head (sort keys %Seq_head_lysogenic){
		if ($head =~ /fragment/){
			delete $Seq_head_lysogenic{$head};
			$Seq_head_prophage{$head} = "prophage";
		}
	}
	
	%Scf2state = (%Seq_head_lytic, %Seq_head_lysogenic, %Seq_head_prophage);
	
	# Phage scaffold length
	my %Scf2length = ();
	open IN, "$img_id/$img_id.a.fna.seq_length.txt";
	while (<IN>){
		chomp;
		if (!/^sequence/){
			my @tmp = split (/\t/);
			my $scf = $tmp[0];
			my $length = $tmp[1];
			$Scf2length{$scf} = $length;
		}
	}
	close IN;
	
	# Store prophage coordinates and scaffold length
	my %Frag2info = (); # $frag => all info
	open IN, "$img_id/VIBRANT_$img_id.a/VIBRANT_results_$img_id.a/VIBRANT_integrated_prophage_coordinates_$img_id.a.tsv";
	while (<IN>){
		chomp;
		if (!/^scaffold/){
			my @tmp = split (/\t/);
			my $frag = $tmp[1];
			my $scf = $tmp[0]; my $scf_length = $Scf2length{$scf};
			shift(@tmp); shift(@tmp);
			push @tmp, $scf_length; # contains 7 elements: protein start | protein stop | protein length | nucleotide start | nucleotide stop | nucleotide length | scaffold length
			$Frag2info{$frag} = join("\t",@tmp);
		}
	}
	close IN;
	
	# Store annotation
	my %Pro2anno = (); # $pro => $anno; here the $pro can contain "fragment" inside
	my @VOG_integrase = qw/VOG00041 VOG15133 VOG20969 VOG02658 VOG04024 VOG01778 VOG02371/;
	my %VOG_integrase = map { $_ => 1 } @VOG_integrase;
	my %Int_pro = (); # Store the integrase protein name, $pro => 1; here the $pro can contain "fragment" inside
	my %Int_scf = (); # Store the integrase scf name, $scf => 1; here the $scf can contain "fragment" inside
	my %Scf2pro = ();  # $scf => $pro separated by "\t" (Both can contain fragment inside)
	open IN, "$img_id/VIBRANT_$img_id.a/VIBRANT_results_$img_id.a/VIBRANT_annotations_$img_id.a.tsv";
	while (<IN>){
		chomp;
		if (!/^pro/){
			my $line = $_;
			my @tmp = split (/\t/);
			my $pro = $tmp[0];
			my $scf = $tmp[1];
			$Pro2anno{$pro} = $line; # This line contains protein	scaffold	KO	AMG	KO name	KO evalue	KO score	KO v-score	Pfam	Pfam name	Pfam evalue	Pfam score	Pfam v-score	VOG	VOG name    VOG evalue	VOG score	VOG v-score
			my $vog = $tmp[13];		 # The order: [0] protein	[1] scaffold	[2] KO	[3] AMG	[4] KO name	[5] KO evalue	[6] KO score	[7] KO v-score	[8] Pfam	[9] Pfam name	[10] Pfam evalue	[11] Pfam score	[12] Pfam v-score	[13] VOG	[14] VOG name    [15] VOG evalue	[16] VOG score	[17] VOG v-score
			if ($vog and exists $VOG_integrase{$vog}){
				$Int_pro{$pro} = 1;
				$Int_scf{$scf} = 1;
			}
			
			
			if (!exists $Scf2pro{$scf}){
				$Scf2pro{$scf} = $pro;
			}else{
				$Scf2pro{$scf} .= "\t".$pro; # Store the protein in a scaffold by order
			}
		}
	}
	close IN;
	
	# Check 2 or over 3 prophages on a scaffold
	my %Prophage_2 = (); # $scf => $frag separated by "\t"; $scf here does not contain fragment inside
	my %Prophage_3 = (); # $scf => $frag separated by "\t"; $scf here does not contain fragment inside
	my %Prophage_4 = (); # $scf => $frag separated by "\t"; $scf here does not contain fragment inside
	my %Scf2frag = (); # $scf => [0] $frag separated by "\t", [1] $frag_num
	foreach my $frag (sort keys %Frag2info){
		my ($scf) = $frag =~ /^(.+?)\_fragment/;  
		if (! exists $Scf2frag{$scf}[0]){
			$Scf2frag{$scf}[0] = $frag;
			$Scf2frag{$scf}[1]++;
		}else{
			$Scf2frag{$scf}[0] .= "\t".$frag;
			$Scf2frag{$scf}[1]++;
		}
	}
	
	foreach my $scf (sort keys %Scf2frag){
		my $frag_num = $Scf2frag{$scf}[1];
		if ($frag_num == 2){
			$Prophage_2{$scf} = $Scf2frag{$scf}[0];
		}elsif ($frag_num == 3){
			$Prophage_3{$scf} = $Scf2frag{$scf}[0];
		}elsif ($frag_num == 4){
			$Prophage_4{$scf} = $Scf2frag{$scf}[0];
		}
	}
	
	my %Combined_phage_scaffold = ();  # $scf => 1
	# For 2 prophages on a scaffold, check to combine or just leave them as original
	foreach my $scf (sort keys %Prophage_2){
		my @Frags = split (/\t/, $Prophage_2{$scf});
		@Frags = sort { $a cmp $b } @Frags; # Sort fragments
		my $frag1 = $Frags[0]; my $frag2 = $Frags[1]; 
		my ($frag1_idx) = $frag1 =~ /^.+\_(\d+?)$/;
		my ($frag2_idx) = $frag2 =~ /^.+\_(\d+?)$/;
		# If two prophages are not adjecent
		if (abs($frag1_idx - $frag2_idx) > 1){
			# Do not do anything
		}else{ # If two prophages are adjecent
			#$Frag2info{$frag} contains 7 elements: protein start | protein stop | protein length | nucleotide start | nucleotide stop | nucleotide length | scaffold length
			my @tmp1 = split (/\t/,$Frag2info{$frag1});
			my $protein_start1 = $tmp1[0];my $protein_stop1 = $tmp1[1];my $protein_length1 = $tmp1[2];
			my $nucleotide_start1 = $tmp1[3];my $nucleotide_stop1 = $tmp1[4];my $nucleotide_length1 = $tmp1[5];
			my $scaffold_length = $tmp1[6];
			
			my @tmp2 = split (/\t/,$Frag2info{$frag2});
			my $protein_start2 = $tmp2[0];my $protein_stop2 = $tmp2[1];my $protein_length2 = $tmp2[2];
			my $nucleotide_start2 = $tmp2[3];my $nucleotide_stop2 = $tmp2[4];my $nucleotide_length2 = $tmp2[5];			
			
			my $protein_distance = $protein_start2 - $protein_stop1 - 1;
			my $nt_distance = $nucleotide_start2 - $nucleotide_stop1 - 1;
			my $front_region = $nucleotide_start1 - 1;
			my $end_region = $scaffold_length - $nucleotide_stop2;
			my $nt_distance_perc = $nt_distance / $scaffold_length;
			my $front_region_perc = $front_region / $scaffold_length;
			my $end_region_perc = $end_region / $scaffold_length;
			my $total_host_region_perc = $nt_distance_perc + $front_region_perc + $end_region_perc;
			if ($protein_distance <= 20 and $nt_distance <= 20000 and $total_host_region_perc <= 0.3){ # Criteria that 2 prophages should be combined
				delete $Scf2state{$frag1}; delete $Scf2state{$frag2};
				
				$Scf2state{$scf} = "lytic"; 
				foreach my $key (sort keys %Int_scf){
					if ($key =~ /$scf\_fragment\_/){
						$Scf2state{$scf} = "lysogenic_but_non-prophage"; 
					}
				}
				
				$Combined_phage_scaffold{$scf} = 1;
			}	
		}
	} # end; for 2 prophages on a scaffold, check to combine or just leave them as original

	# for 3 prophages on a scaffold, check to combine or just leave them as original
	foreach my $scf (sort keys %Prophage_3){ # Right bracket for 3 prophages on a scaffold 
		my @Frags = split (/\t/, $Prophage_3{$scf});
		@Frags = sort { $a cmp $b } @Frags; # Sort fragments
		my $frag1 = $Frags[0]; my $frag2 = $Frags[1]; my $frag3 = $Frags[2]; 
		
		my @tmp1 = split (/\t/,$Frag2info{$frag1});
		my $protein_start1 = $tmp1[0];my $protein_stop1 = $tmp1[1];my $protein_length1 = $tmp1[2];
		my $nucleotide_start1 = $tmp1[3];my $nucleotide_stop1 = $tmp1[4];my $nucleotide_length1 = $tmp1[5];
		my $scaffold_length = $tmp1[6];
		
		my @tmp2 = split (/\t/,$Frag2info{$frag2});
		my $protein_start2 = $tmp2[0];my $protein_stop2 = $tmp2[1];my $protein_length2 = $tmp2[2];
		my $nucleotide_start2 = $tmp2[3];my $nucleotide_stop2 = $tmp2[4];my $nucleotide_length2 = $tmp2[5];	

		my @tmp3 = split (/\t/,$Frag2info{$frag3});
		my $protein_start3 = $tmp3[0];my $protein_stop3 = $tmp3[1];my $protein_length3 = $tmp3[2];
		my $nucleotide_start3 = $tmp3[3];my $nucleotide_stop3 = $tmp3[4];my $nucleotide_length3 = $tmp3[5];	
	
		my $protein_distance_1_2 = $protein_start2 - $protein_stop1 - 1;
		my $nt_distance_1_2 = $nucleotide_start2 - $nucleotide_stop1 - 1;
		my $protein_distance_2_3 = $protein_start3 - $protein_stop2 - 1;
		my $nt_distance_2_3 = $nucleotide_start3 - $nucleotide_stop2 - 1;
		my $front_region = $nucleotide_start1 - 1;
		my $end_region = $scaffold_length - $nucleotide_stop3;
		
		my $nt_distance_1_2_perc = $nt_distance_1_2 / $scaffold_length;
		my $nt_distance_2_3_perc = $nt_distance_2_3 / $scaffold_length;
		my $front_region_perc = $front_region / $scaffold_length;
		my $end_region_perc = $end_region / $scaffold_length;
		my $total_host_region_perc = $nt_distance_1_2_perc + $nt_distance_2_3_perc + $front_region_perc + $end_region_perc;	

		if ($protein_distance_1_2 <= 20 and $protein_distance_2_3 <= 20 and $nt_distance_1_2 <= 20000 and $nt_distance_2_3 <= 20000 and $total_host_region_perc <= 0.3){ # Criteria that 3 prophages should be combined
				delete $Scf2state{$frag1}; delete $Scf2state{$frag2}; delete $Scf2state{$frag3};
				
				$Scf2state{$scf} = "lytic"; 
				foreach my $key (sort keys %Int_scf){
					if ($key =~ /$scf\_fragment\_/){
						$Scf2state{$scf} = "lysogenic_but_non-prophage"; 
					}
				}	

				$Combined_phage_scaffold{$scf} = 1;				
		}	
	}# Left bracket for 3 prophages on a scaffold 
	
	# For 4 prophages on a scaffold, check to combine or just leave them as original
	foreach my $scf (sort keys %Prophage_4){ # Right bracket for 4 prophages on a scaffold 
		my @Frags = split (/\t/, $Prophage_4{$scf});
		@Frags = sort { $a cmp $b } @Frags; # Sort fragments
		my $frag1 = $Frags[0]; my $frag2 = $Frags[1]; my $frag3 = $Frags[2]; my $frag4 = $Frags[3]; 
		
		my @tmp1 = split (/\t/,$Frag2info{$frag1});
		my $protein_start1 = $tmp1[0];my $protein_stop1 = $tmp1[1];my $protein_length1 = $tmp1[2];
		my $nucleotide_start1 = $tmp1[3];my $nucleotide_stop1 = $tmp1[4];my $nucleotide_length1 = $tmp1[5];
		my $scaffold_length = $tmp1[6];
		
		my @tmp2 = split (/\t/,$Frag2info{$frag2});
		my $protein_start2 = $tmp2[0];my $protein_stop2 = $tmp2[1];my $protein_length2 = $tmp2[2];
		my $nucleotide_start2 = $tmp2[3];my $nucleotide_stop2 = $tmp2[4];my $nucleotide_length2 = $tmp2[5];	

		my @tmp3 = split (/\t/,$Frag2info{$frag3});
		my $protein_start3 = $tmp3[0];my $protein_stop3 = $tmp3[1];my $protein_length3 = $tmp3[2];
		my $nucleotide_start3 = $tmp3[3];my $nucleotide_stop3 = $tmp3[4];my $nucleotide_length3 = $tmp3[5];	
		
		my @tmp4 = split (/\t/,$Frag2info{$frag4});
		my $protein_start4 = $tmp4[0];my $protein_stop4 = $tmp4[1];my $protein_length4 = $tmp4[2];
		my $nucleotide_start4 = $tmp4[3];my $nucleotide_stop4 = $tmp4[4];my $nucleotide_length4 = $tmp4[5];			
	
		my $protein_distance_1_2 = $protein_start2 - $protein_stop1 - 1;
		my $nt_distance_1_2 = $nucleotide_start2 - $nucleotide_stop1 - 1;
		my $protein_distance_2_3 = $protein_start3 - $protein_stop2 - 1;
		my $nt_distance_2_3 = $nucleotide_start3 - $nucleotide_stop2 - 1;
		my $protein_distance_3_4 = $protein_start4 - $protein_stop3 - 1;
		my $nt_distance_3_4 = $nucleotide_start4 - $nucleotide_stop3 - 1;		
		
		my $front_region = $nucleotide_start1 - 1;
		my $end_region = $scaffold_length - $nucleotide_stop4;
		
		my $nt_distance_1_2_perc = $nt_distance_1_2 / $scaffold_length;
		my $nt_distance_2_3_perc = $nt_distance_2_3 / $scaffold_length;
		my $nt_distance_3_4_perc = $nt_distance_3_4 / $scaffold_length;
		my $front_region_perc = $front_region / $scaffold_length;
		my $end_region_perc = $end_region / $scaffold_length;
		my $total_host_region_perc = $nt_distance_1_2_perc + $nt_distance_2_3_perc + $nt_distance_3_4_perc + $front_region_perc + $end_region_perc;	

		if ($protein_distance_1_2 <= 20 and $protein_distance_2_3 <= 20 and $protein_distance_3_4 <= 20 and $nt_distance_1_2 <= 20000 and $nt_distance_2_3 <= 20000 and $nt_distance_3_4 <= 20000 and $total_host_region_perc <= 0.3){ # Criteria that 4 prophages should be combined
				delete $Scf2state{$frag1}; delete $Scf2state{$frag2}; delete $Scf2state{$frag3}; delete $Scf2state{$frag4};
				
				$Scf2state{$scf} = "lytic"; 
				foreach my $key (sort keys %Int_scf){
					if ($key =~ /$scf\_fragment\_/){
						$Scf2state{$scf} = "lysogenic_but_non-prophage"; 
					}
				}	
				
				$Combined_phage_scaffold{$scf} = 1;
		}	
	}# Left bracket for 4 prophages on a scaffold 

	# Print phage result tsv and grab phage sequences
	# Only keep phage results that pass min length
	`mkdir $img_id/VIBRANT_$img_id.a.v2.min$min_length`;
	my $result_dir = "$img_id/VIBRANT_$img_id.a.v2.min$min_length";
	my %Seq_fna = _store_seq("$img_id/$img_id.a.fna");
	my %Seq_fna_former_combined = _store_seq("$img_id/VIBRANT_$img_id.a/VIBRANT_phages_$img_id.a/$img_id.a.phages_combined.fna");
	my %Seq_faa = _store_seq("$img_id/VIBRANT_$img_id.a/$img_id.a.prodigal.faa");
	my %Seq_faa_former_combined = _store_seq("$img_id/VIBRANT_$img_id.a/VIBRANT_phages_$img_id.a/$img_id.a.phages_combined.faa");
	my %Seq_ffn = _store_seq("$img_id/VIBRANT_$img_id.a/$img_id.a.prodigal.ffn");
	my %Seq_ffn_former_combined = _store_seq("$img_id/VIBRANT_$img_id.a/VIBRANT_phages_$img_id.a/$img_id.a.phages_combined.ffn");
	
	# Print out phage scaffolds
	my %Seq_fna_out = (); my %Seq_faa_out = (); my %Seq_ffn_out = ();
	%Seq_fna_out = %Seq_fna_former_combined;
	%Seq_faa_out = %Seq_faa_former_combined;
	%Seq_ffn_out = %Seq_ffn_former_combined;
	
	# Change output fna hash
	foreach my $key (sort keys %Seq_fna_out){ # Delete old fragments
		foreach my $scf (sort keys %Combined_phage_scaffold){
			if ($key =~ /$scf\_fragment/){
				delete $Seq_fna_out{$key};
			}
		}
	}
	
	foreach my $scf (sort keys %Combined_phage_scaffold){ # Add new whole scaffold
		$Seq_fna_out{">$scf"} = $Seq_fna{">$scf"};
	}
	
	# Change output faa hash
	foreach my $key (sort keys %Seq_faa_out){ # Delete old fragment proteins
		foreach my $scf (sort keys %Combined_phage_scaffold){
			if ($key =~ /$scf\_fragment/){
				delete $Seq_faa_out{$key};
			}
		}
	}
	
	foreach my $scf (sort keys %Combined_phage_scaffold){ # Add new whole scaffold
		foreach my $key (sort keys %Seq_faa){
			if ($key =~ /$scf\_/){
				$Seq_faa_out{$key} = $Seq_faa{$key};
			}
		}
	}
	
	# Change output ffn hash
	foreach my $key (sort keys %Seq_ffn_out){ # Delete old fragment proteins
		foreach my $scf (sort keys %Combined_phage_scaffold){
			if ($key =~ /$scf\_fragment/){
				delete $Seq_ffn_out{$key};
			}
		}
	}
	
	foreach my $scf (sort keys %Combined_phage_scaffold){ # Add new whole scaffold
		foreach my $key (sort keys %Seq_ffn){
			if ($key =~ /$scf\_/){
				$Seq_ffn_out{$key} = $Seq_ffn{$key};
			}
		}
	}	
	
	# Delete sequences with less than $min_length in hash of fna
	foreach my $key (sort keys %Seq_fna_out){
		my $len = length ($Seq_fna_out{$key});
		if ($len < $min_length){
			delete $Seq_fna_out{$key};
		}
	}
	
	# Print phage and prophage fna
	open OUT, ">$result_dir/$img_id.a.phage_scaffold.fna";
	open OUT2, ">$result_dir/$img_id.a.prophage_scaffold.fna";
	foreach my $key (sort keys %Seq_fna_out){
		print OUT "$key\n$Seq_fna_out{$key}\n";
		if ($key =~ /fragment/){
			print OUT2 "$key\n$Seq_fna_out{$key}\n";
		}
	}
	close OUT;
	close OUT2;

	# Print phage and prophage faa
	open OUT, ">$result_dir/$img_id.a.phage_scaffold.faa";
	open OUT2, ">$result_dir/$img_id.a.prophage_scaffold.faa";
	foreach my $key (sort keys %Seq_faa_out){
		my ($scf) = $key =~ /^>(.+)\_(\d+?)$/;
		if (exists $Seq_fna_out{">$scf"}){ # Only write faa that are from > min_length scaffolds
			print OUT "$key\n$Seq_faa_out{$key}\n";
			if ($key =~ /fragment/){
				print OUT2 "$key\n$Seq_faa_out{$key}\n";
			}		
		}
	}
	close OUT;
	close OUT2;
	
	# Print phage and prophage ffn
	open OUT, ">$result_dir/$img_id.a.phage_scaffold.ffn";
	open OUT2, ">$result_dir/$img_id.a.prophage_scaffold.ffn";
	foreach my $key (sort keys %Seq_ffn_out){
		my ($scf) = $key =~ /^>(.+)\_(\d+?)$/;
		if (exists $Seq_fna_out{">$scf"}){ # Only write ffn that are from > min_length scaffolds
			print OUT "$key\n$Seq_ffn_out{$key}\n";
			if ($key =~ /fragment/){
				print OUT2 "$key\n$Seq_ffn_out{$key}\n";
			}		
		}
	}
	close OUT;
	close OUT2;	
	
	# Print phage result
	open OUT, ">$result_dir/phage_results.min$min_length.txt";
	foreach my $scf (sort keys %Scf2state){
		if ($Seq_fna_out{">$scf"}){
			my $length_scf = length ($Seq_fna_out{">$scf"});
			print OUT "$scf\t$Scf2state{$scf}\t$length_scf\n";
		}
	}
	close OUT;
	
	# Print prophage coordinates
	open OUT, ">$result_dir/prophage_coordinates.txt";
	print OUT "scaffold\tfragment\tprotein start\tprotein stop\tprotein length\tnucleotide start\tnucleotide stop\tnucleotide length\n";
	foreach my $scf (sort keys %Scf2state){
		if ($scf =~ /fragment/){
			my $length_scf = length ($Seq_fna_out{">$scf"});
			if ($length_scf >= $min_length){
				my @tmp = split (/\t/, $Frag2info{$scf});
				pop @tmp; # remove the last one element;
				# scaffold        fragment        protein start   protein stop    protein length  nucleotide start        nucleotide stop nucleotide length
				my $fragment = $scf;
				my ($scaffold) = $fragment =~ /^(.+?)\_fragment/;
				my $line = join ("\t",@tmp);
				print OUT "$scaffold\t$fragment\t$line\n";
			}
		}
	}
	close OUT;
}



# Subroutine

sub _store_seq_head_to_state{
	my $file = $_[0];
	my $state = $_[1];
	my %Seq_head = ();
	open _IN, "$file";
	while (<_IN>){
		chomp;
		if (/^>/){
			my $head = $_;
			$head =~ s/^>//g;
			$Seq_head{$head} = $state;
		}
	}
	close _IN;
	return %Seq_head;
}	

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