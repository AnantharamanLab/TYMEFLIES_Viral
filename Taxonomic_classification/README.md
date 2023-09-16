# Taxonomic classification

(The method was based on *Nucleic Acids Res*. 2021 Jan 8;49(D1):D764-D775. doi: 10.1093/nar/gkaa946. Link: https://academic.oup.com/nar/article/49/D1/D764/5952208)

**1 Run diamond to NCBI RefSeq to determine phage genome classification**

Use all phage faa files to compare against NCBI Viral RefSeq proteins with diamond.

To see whether >= 30% of the proteins for a faa (phage genome) have a hit to Viral RefSeq, only count phage genomes that are above this cutoff.

Get the consensus taxonomic affiliation based on the best hits of individual proteins (>= 50% majority rule) for a phage genome.

[script] 05.Taxonomic_classification.step1.run_diamond_to_NCBI_RefSeq_viral.pl

**2 Run hmmsearch to VOG marker HMM to determine phage genome classification**

Use all phage faa files to compare against VOG 587 marker HMM with hmmsearch (Note VOG97 is used here, since the VOG marker IDs are corresponding to VOG IDs in VOG97, later version VOG database might cause incongruencies on VOG IDs).

Get the consensus taxonomic affiliation based on individual markers detected by a simple plurality rule, If there are multiple conflicting markers detected .

[script] 05.Taxonomic_classification.step2.run_hmmsearch_to_VOG_marker.pl

**3** **Use geNomad to annotate all viral scaffolds to get taxonomy**

First we used 1000 Ns to link multiple fasta sequences to make multiple-fasta viruses to be single-fasta viruses (to meet the input requirement of geNomad). Then we used the default settings of "geNomad annotate" to annotate viruses and get taxonomy from the annotation result.

**4 Integrate all results and provide the final taxonomic classification result**

Store the taxonomic classification results from the above three steps ("NCBI RefSeq viral protein searching", "VOG marker HMM searching", and "geNomad classifying"). And then get taxonomy based on other members' taxonomy (only the hits of "NCBI RefSeq viral protein searching" will be counted) from each genus and expand the tax to all the members within this genus (Genus-level vOTU LCA assigning). Last, write down the final taxonomic classification result.

[script] 05.Taxonomic_classification.step4.integrate_all_results.pl





# Script for wide usage

Since that these two scripts are very custom for input datasets in this study, we modified these scripts for wide usage. 

**Usage:** The input should be two files:
        1) all virus protein sequences (file ended with "faa")
        2) a map file of protein to genome; this map file should contain two columns: 1. each protein 2. the corresponding genome of this protein
        An example here:
         protein_id	genome_id
		 protein_A	genome_A
		 protein_B	genome_A
		 protein_C	genome_A
		 protein_D	genome_B
		 protein_E	genome_B
		 protein_F	genome_B	

perl 05.Taxonomic_classification.step1.run_diamond_to_NCBI_RefSeq_viral.pl [virus_protein_sequence.faa] [map_file.txt] 
The first two items are the inputs, the 3rd item is the output taxonomy result

perl 05.Taxonomic_classification.step2.run_hmmsearch_to_VOG_marker.pl [virus_protein_sequence.faa] [map_file.txt]
The first two items are the inputs, the 3rd item is the output taxonomy result

The perl script will generate a folder called "Taxonomic_classification" in the current directory.



**Auxiliary script:**

[script] pre01.generate_protein_to_genome_map_for_single_contig_virus.pl

This script will help you to generate a protein to genome map file for single-contig virus.

Usage: The input should be the all virus protein file (ended with "faa")
perl pre01.generate_protein_to_genome_map_for_single_contig_virus.pl [virus_protein_sequence.faa] [map_file.txt] 
The first item is the input, the second item is the output (the map file)

Note: The input faa file should be a prodigal-formate protein sequence file
Note: Note that the virus genome is made up of just one contig/scaffold



