# Time-series analysis

(Include AMG coverage ratio analysis and micro-/macrodiversity analysis for viruses using time-series metagenomes across 20 years)

**1 Get the information of high occurrence species (>= 40)**

The occurrence info is based on the input file of "Species2occurrence_n_abundance.txt" (generated from Step 3 in [Mapping metagenomic assemblies](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Mapping_metagenomic_assemblies)). This occurrence is based on directly binned viral genomes in individual metagenomes. 

The information contains:1) Tax 2) Host tax 3) AMG KO within. This script also contains a filtering step: Only keep species that contain >= 10 genomes.

[script] 08.Time_series_analysis.step1.get_information_of_high_occurrence_species.pl

**2 Map all metagenomic reads to the collection of species representatives**

Map all metagenomic reads to the collection of the representative genomes from individual species. We only included representative genomes that carry at least one AMG.

[script] 08.Time_series_analysis.step2.map_metagenomic_reads_to_the_collection_of_species_representatives.pl

**3 Grep gene files (in prodigal format) for "All_phage_species_rep_gn_containing_AMG.fasta"**

"All_phage_species_rep_gn_containing_AMG.fasta" is the fasta file contains species representative genomes (with at least one AMG) generated from the last step (Step 2). We parsed the ffn file (the prodigal-annotation result containing all genes) from each VIBRANT result folder to get the new gene headers and the corresponding sequences. 

[script] 08.Time_series_analysis.step3.grep_gene_files_for_all_rep_gn_containing_AMG.pl

**4 Run MetaPop**

1) "--id_min" set to 92 (3% lower than the exact species genome boundary as suggested)
2) Use genes provided by me (generated from Step 3)
3) "--genome_detection_cutoff" set to 70 (Percent of bases that must be covered for a sequence to be considered detected for macrodiversity analyses - 70%)

4) "--snp_scale" set to "both" (Including both SNPs detected for a genome in each sample alone, and SNPs for that genome across all samples.)

[script] 08.Time_series_analysis.step4.run_metapop.pl

**5 Parse MetaPop result to get AMG coverage**

In this script, two AMG coordinate files need to be provided:

1) For non-prophage viruses: AMG coordinate file was parsed from the gene file (generated from Step 3).

2) For prophage: AMG coordinate file was the concatenated file of all "VIBRANT_integrated_prophage_coordinates_*.a.tsv" from VIBRANT result folders.

A custom Python 3 script "cov_by_region.py" was used to get AMG coverage. Note that this script should be run under conda env "python_scripts_env_Jan2022.yml". Both the script and yml file were provided here.

The script uses the outputs from Directory "04.Depth_per_Pos" generated from MetaPop.

In the resulted file "*.AMG_cov.txt", five columns were included:

1) scaffold (scaffold ID)

2) region (gene ID of AMG)

3) avg whole scaffold #1 (mean coverage of scaffold)

4) avg partial scaffold #1 (mean coverage of scaffold excluding AMG regions)

5) avg region #1 (mean coverage of AMG)

[script] 08.Time_series_analysis.step5.parse_metapop_result_to_get_AMG_coverage.pl and cov_by_region.py

**6 Parse the AMG coverage**

Processing the coverage result with the following criteria:

1) screen bin with its any scaffold with < 0.01 coverage 

2) normalize the coverage by read numbers (100M reads per metagenome)

We have processed both the viral genome coverages and AMG coverage ratios (AMG coverage / viral genome coverage).

[script] 08.Time_series_analysis.step6.parse_AMG_coverage.pl

**7 Get viral genome and AMG variation**

Parse to get the viral genome coverage and AMG coverage ratio variation pattern across all metagenomes. Process the viral genome and assign it as "present" in a metagenome with both coverage >= 1.5 and breadth >= 70%.

The resulted file "AMG_gene_cov_ratio_variation_table.txt" contains the following columns:

1) KO

2) KO detail

3) Viral genome

4) Distribution of viral genome (present in how many metagenomes)

5) Coverage mean across all metagenomes

6) Distribution of AMG (present in how many metagenomes)

7) - 471) AMG coverage ratio in 475 metagenomes

[script] 08.Time_series_analysis.step7.get_viral_gn_and_AMG_variation.pl

**8 Get AMG coverage ratio distribution in each month and year-month**

Note that we only included AMGs with distribution >= 233 (0.5 * 465, that is half of the metagenome number). In this script, we calculated both AMG coverage ratio distribution in each month and year-month and AMG-containing viral genome coverage distribution in each month and year-month.

[script] 08.Time_series_analysis.step8.get_AMG_cov_ratio_variation_in_month_and_year_month.pl

**9 Visualize AMG coverage ratio and corresponding viral genome coverage variations**

The script generates two figures:

1) AMG coverage ratio figure

It contains three subpanels of A) heatmap of the number of metagenomes for each month, B) 20 facets of AMG coverage ratio plots for each year-month, and C) AMG coverage ratio plot for each month (aggregated across 20 years).

2) AMG-containing viral genome coverage figure

It contains three subpanels of A) heatmap of the number of metagenomes for each month, B) 20 facets of AMG-containing viral genome coverage plots for each year-month, and C) AMG-containing viral genome coverage plot for each month (aggregated across 20 years).

08.Time_series_analysis.step9.visualize_AMG_and_corresponding_viral_gn_cov_ratio_varition.R