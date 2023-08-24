#!/usr/bin/perl

use strict;
use warnings;

# AIM: Move all TYMEFLIES_MAGs that passed CheckM completeness and contamination requirements and use Prodigal to get their faa files


# Step 1 Move all MAGs
my $work_dir_TYMEFLIES_MAGs = "/storage1/data11/TYMEFLIES_phage/Robin_MAGs";
=pod
`mkdir $work_dir_TYMEFLIES_MAGs/All_passed_MAGs`;
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
			`cp $MAG_addr $work_dir_TYMEFLIES_MAGs/All_passed_MAGs`;
		}
	}
}
close IN;
=cut

# Step 2 Annotate all MAGs using Prodigal
my $input_genome_folder = $work_dir_TYMEFLIES_MAGs."/All_passed_MAGs";

## Find all FASTA files in the input genome folder
my @fasta_files = glob("$input_genome_folder/*.fasta");

# Open the output file for writing
open(my $output_fh, '>', 'tmp.run_prodigal.sh') or die "Cannot open output file: $!";
foreach my $fasta_file (@fasta_files) {    
    my ($gn_id) = $fasta_file =~ /\/([^\/]+)\.fasta$/; # Extract the file name without the extension    
    my $prodigal_cmd = "prodigal -i $fasta_file -a $input_genome_folder/$gn_id.faa -o $input_genome_folder/$gn_id.gff -f gff -p meta -q"; # Define the Prodigal command    
    print $output_fh "$prodigal_cmd\n"; # Print the Prodigal command to the output file
}
close($output_fh);
=pod
# Run the script using parallel
system('cat tmp.run_prodigal.sh | parallel -j 20') == 0
    or die "Error running the script with parallel: $!";
=cut


