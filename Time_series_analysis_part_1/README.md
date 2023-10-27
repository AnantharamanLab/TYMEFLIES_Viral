# Time-series analysis - Part 1 AMG ratio and viral genome coverage analysis

(Include AMG ratio and viral genome coverage analysis using time-series metagenomes across 20 years)

**1 Get the information of high occurrence species (>= 40)**

The occurrence info is based on the input file of "Species2occurrence_n_abundance.txt" (generated from Step 3 in [Mapping metagenomic assemblies](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Mapping_metagenomic_assemblies)). This occurrence is based on directly binned viral genomes in individual metagenomes. 

The information contains:1) Tax 2) Host tax 3) AMG KO within. This script also contains a filtering step: Only keep species that contain >= 10 genomes.

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step1.get_information_of_high_occurrence_species.pl

**2 Get the AMG counterpart genes and flankings from all metagenomes**

It should include AMG counterpart genes and their flanking regions (+-150bp) in the mapping reference for the next step. We first run hmmsearch of all metagenome proteins to AMG KOs, and then grab the AMG counterpart genes and their flankings from positive hmmsearch hits. We excluded all viral scaffolds in the final AMG counterpart gene located scaffolds.

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step2.get_AMG_counterpart_gene_and_flankings.pl

**3 Map all metagenomic reads to the collection of species representatives and AMG counterpart genes and flankings**

Map all metagenomic reads to the collection of species representatives and AMG counterpart genes and flankings.  The resulting sam files were converted to bam files and subjected to read identity filtering with a cutoff threshold of 90%.

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step3.map_metagenomic_reads_to_the_collection_of_species_representatives.pl

**4 Grep gene files (in prodigal format) for "All_phage_species_rep_gn.fasta"**

"All_phage_species_rep_gn.fasta" is the fasta file contains all the species representative genomes. We parsed the ffn file (the prodigal-annotation result containing all genes) from each VIBRANT result folder to get the new gene headers and the corresponding sequences. 

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step4.grep_gene_files_for_all_rep_gn_containing_AMG.pl

**5 Filter bam files using only the collection of viral species representative genomes as the reference**

All the "*.viral_species_rep.id90.bam" files from the Step 3 were placed into a new folder. Filter bam files by scaffold names of viral species representative genomes. 

A custom Python 3 script "filter_bam_by_reference.py" was used to filter bam. Note that this script should be run under conda env "python_scripts_env_Jan2022.yml" (requirement: pysam >= 0.16). Both the script and yml file were provided here.

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step5.filter_bam_files.pl

**6 Run MetaPop**

1) "--id_min" set to 93 (species genome boundary, this value was suggested here  - the section of "Reducing read mis-mapping by adjusting min_read_ani" in https://instrain.readthedocs.io/en/latest/important_concepts.html#strain-level-comparisons-and-popani)
2) Use genes provided by me (generated from Step 3)
3) "--genome_detection_cutoff" set to 70 (Percent of bases that must be covered for a sequence to be considered detected for macrodiversity analyses - 70%)

4) "--snp_scale" set to "both" (Including both SNPs detected for a genome in each sample alone, and SNPs for that genome across all samples.)

Note:  The custom gene file should be strictly in the same format/shape as a Prodigal-generated one. The gene ID should be re-ordered (naturally ordered) if you want to use a custom gene file. For instance, 11 is in front of 2 , 3, and 4, which needs to be re-ordered. Meanwhile, the gene ID should start from 1, and the start and stop points of each gene in a scaffold should start from the very beginning of the scaffold (but not the intermediate place of a scaffold). A custom perl script (Side23.change_genes_headers.pl) was provided to change the headers of the gene file.

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step6.run_metapop.pl

**7 Parse MetaPop result to get AMG coverage**

In this script, two AMG coordinate files need to be provided:

1) For non-prophage viruses: AMG coordinate file was parsed from the gene file (generated from Step 4).

2) For prophage: AMG coordinate file was the concatenated file of all "VIBRANT_integrated_prophage_coordinates_*.a.tsv" from VIBRANT result folders.

A custom Python 3 script "cov_by_region.py" was used to get AMG coverage. Note that this script should be run under conda env "python_scripts_env_Jan2022.yml". Both the script and yml file were provided here.

The script uses the outputs from Directory "04.Depth_per_Pos" generated from MetaPop.

In the resulted file "*.AMG_cov.txt", five columns were included:

1) scaffold (scaffold ID)

2) region (gene ID of AMG)

3) avg whole scaffold #1 (mean coverage of scaffold)

4) avg partial scaffold #1 (mean coverage of scaffold excluding AMG regions)

5) avg region #1 (mean coverage of AMG)

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step7.parse_metapop_result_to_get_AMG_coverage.pl

**8 Parse the AMG coverage**

Processing the coverage result with the following criteria:

1) screen bin with its any scaffold with < 0.01 coverage 

2) normalize the coverage by read numbers (100M reads per metagenome)

We have processed both the viral genome coverages and AMG coverage ratios (AMG coverage / viral genome coverage).

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step8.parse_AMG_coverage.pl

**9 Get viral genome and AMG variation**

Parse to get the viral genome coverage and AMG coverage ratio variation pattern across all metagenomes. Process the viral genome and assign it as "present" in a metagenome with both coverage >= 0.33 and breadth >= 50%.

The resulted file "AMG_gene_cov_ratio_variation_table.txt" contains the following columns:

1. Col 1: KO

2. Col 2: KO detail

3. Col 3: Viral genome

4. Col 4: Distribution of viral genome (present in how many metagenomes)

5. Col 5: Coverage mean across all metagenomes

6. Col 6: Distribution of AMG (present in how many metagenomes)
7. Col 7-476: AMG coverage ratio in 471 metagenomes

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step9.get_viral_gn_and_AMG_variation.pl

**10 Get AMG coverage ratio distribution in each season and year-season**

Note that we only included AMGs with distribution >= 5. In this script, we calculated both AMG coverage ratio distribution in each season and year-season and AMG-containing viral genome coverage distribution in each season and year-season.

"MetaPop/KO2season_ko_details.txt" is a file that contains two columns: "KO ID" and "KO details"

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step10.get_AMG_cov_ratio_variation_in_season_and_year_season.pl

**14 Map all metagenomic reads from each year to the collection of species representatives and AMG counterpart genes and flankings**

Map all metagenomic reads from each year to the collection of the representative genomes from individual species and AMG counterpart genes and flankings. 

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step14.map_metagenomic_reads_from_each_year_to_the_collection_of_species_representatives.pl

**15 Filter bam files using only the collection of viral species representative genomes containing AMGs as the reference**

All the "*.viral_species_rep.id90.bam" files from the Step 14 were placed into a new folder. Filter bam files by scaffold names of viral species representative genomes containing AMGs. 

A custom Python 3 script "filter_bam_by_reference.py" was used to filter bam. Note that this script should be run under conda env "python_scripts_env_Jan2022.yml". Both the script and yml file were provided here.

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step15.filter_bam_files.for_each_year.pl

**16 Run MetaPop for bam files from each year**

The settings were the same as those listed in Step 6. The mapping reference is "All_phage_species_rep_gn_containing_AMG.fasta".

[script] 

08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step16.run_metapop.for_each_year.pl

**17 Calculate the AMG containing viral genome statistics**

This script is written to calculate the following statistics of AMG containing viral genomes:

(1) The percentage of species representatives containing AMG 

(2) The psbA-containing species rep to species with any members containing psbA ratio

(3) The virus completeness to specific AMG-containing species members AMG containing percentage

​     To see if the virus completeness degree will influence the AMG-containing situation of the species members. The following AMGs are inspected: 

```sh
# AMG gene: KO ID
'psbA': 'K02703',
'psbD': 'K02706',
'pmoC': 'K10946',
'ahbD': 'K22227',
'katG': 'K03782',
'gpx': 'K00432',
'cobS': 'K09882',
'cobT': 'K09883',
'nadE': 'K01916',
'cysC': 'K00860',
'cysD': 'K00957',
'cysH': 'K00390',
'cysK': 'K01738'
```
[script]

08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step17.AMG_containing_gn_stat.py

**18 Visualize the virus completeness to specific AMG-containing species members AMG containing percentage**

For each AMG, we plotted two bars representing two genome completeness ranges: 75-100% and 50-75%. Each bar represents the percentage of AMG containing percentage of the species members that fall into the corresponding genome completeness range.

[script]

08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step18.AMG_containing_gn_stat.visualize.R

**19 Parse to get other viral genome (no-AMG containing viral genome) abundance**

Parse to get other viral genome abundance (viruses that do not contain any AMGs)
Screen viral genome with the following two criteria:
(1) Screen viral genome with its any scaffold with < 0.01 coverage 
(2) Process the viral genome presence with both coverage (>= 0.33) and breadth (>= 50%)

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step19.parse_to_get_other_viral_gn_abundance.py

**20 Make the composition pattern table based on viral genome abundance at the family level**

The composition pattern table comprised of the abundances of viral family across all the samples. The sample to season information was also included.

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step20.make_composition_pattern_table.py

**21 Visualize the composition pattern**

Visualize the composition pattern of viral community at the family level across six seasons by NMDS plots and conduct ANOSIM tests for all six seasons. At the same time, visualize the composition pattern of viral community of the seasons of (1) "Fall", "Ice-on", and "Spring", (2) "Spring", "Clearwater", and "Early Summer", (3) "Early Summer", "Late Summer", and "Fall", separately; and conduct ANOSIM tests accordingly.

[script] 08.Time_series_analysis.part1.AMG_ratio_and_viral_gn_analysis.step21.visualize_composition_pattern.R




