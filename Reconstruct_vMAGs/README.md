# Reconstruct vMAGs

**1 Run vRhyme to reconstruct vMAGs from phage scaffolds identified by VIBRANT**

Use default settings. vRhyme link: https://github.com/AnantharamanLab/vRhyme

[script] 03.Reconstruct_vMAGs.step1.run_vRhyme.pl

**2 Find the best bin set from vRhyme result**

Four criteria were used to refine the best bin set from vRhyme result:                     

1) Prophage identified by VIBRANT should be excluded from binning                        

2) Two or more lysogenic (non-prophage) phage scaffolds can not be in the same bin                          

3) Phage scaffolds identified by CheckV as "Complete" should be excluded from binning    

4) The maximum number of bin redundancy should be <= 1  

The resulted faa, ffn, fasta files for each vRhyme bin and vRhyme unbinned genome were provided

[script] 03.Reconstruct_vMAGs.step2.make_new_vRhyme_bins.pl

**3 Run CheckV for each phage genome**

Due to CheckV only allows for single-contig phage genome, we firstly link bin scaffolds with multiple Ns (1500 Ns as suggested by default) to make temporary single-contig phage genome. Then we run CheckV for temporary single-contig phage genome with default settings.

[script] 03.Reconstruct_vMAGs.step3.run_checkV_for_each_phage_genome.pl

[script] path/to/folder/vRhyme_v1.0.0/vRhyme/aux/link_bin_sequences.py (for linking bin scaffolds)

This script can be found in "https://github.com/AnantharamanLab/vRhyme/tree/master/vRhyme/aux"

**4 Summary CheckV results for each phage genome**

[script] 13.Reconstruct_vMAGs.step4.summarize_checkV_result_for_each_phage_genome.pl

**5 & 6 Summarize all AMGs from individual metagenomes** 

Step 5: Generate AMG table for all metagenomes. Use the AMG annotated by VIBRANT and change the protein according to the IMG ID and phage genome ID

Step 6: Add metagenome info and KEGG pathway info to the AMG table. The following information has been added:

Metagenome date and season ('date | date in the year | season')

Metabolisms (mutiple metabolisms separated by '|')

Pathways (mutiple pathways separated by '|')

This script will use the file ["VIBRANT_KEGG_pathways_summary.tsv"](https://github.com/AnantharamanLab/VIBRANT/blob/master/files/VIBRANT_KEGG_pathways_summary.tsv) from VIBRANT.  

[script] 03.Reconstruct_vMAGs.step5.generate_AMG_summary_table.pl and 03.Reconstruct_vMAGs.step6.add_metagenome_and_KEGG_pathway_info_to_AMG_summary_table.pl

**7 Summarize AMG frequencies and abundance monthly and calculate AMG abundance for each year-month (for example, "2000-01")**

This script will generate following results:

1) KO2month2freq.txt

AMG frequencies (normalized by read number per metagenome, 100M reads per metagenome, and metagenome number per month) for each month

2) KO2month2abun.txt

AMG abundances (nnormalized by read number per metagenome, 100M reads per metagenome, and metagenome number per month) for each month

3) KO2year_month2abun.txt

AMG abundances (normalized by read number per metagenome, 100M reads per metagenome, and metagenome number per year-month) for each year-month

4) KO2month_ko_details.txt

Two columns containing KO ID and KO detail for all KOs that appear in all metagenomes

5) Month2num_of_metagenomes.txt

the number of metagenomes in each month 

6) Year_month2num_of_metagenomes.txt

the number of metagenomes in each year_month 

[script] 03.Reconstruct_vMAGs.step7.summarize_AMG_frequency_and_abundance.pl

**8 Grep all AMG proteins and run against dbCAN2 and MEROPS databases to get annotated**

Get all AMG proteins to annotate by dbCAN2 and MEROPS databases.

This will add more columns to "AMG_summary.txt":

"CAZyme ID"

"Enzyme (CAZyme)"

"Substrate (CAZyme)"

"Class of Substrate (CAZyme)"

"MEROPS ID"

"MEROPS details"

[script] 03.Reconstruct_vMAGs.step8.annnoate_AMG_by_dbCAN2_and_MEROPS.pl

**9 Make AMG trend R plots for each month and each year-month**

Use the results generated from Step 6 ("KO2month2abun.txt", "KO2year_month2abun.txt", and "KO2month_ko_details.txt") to draw AMG trend R plots (AMG KO abundance) for each month (aggregated across 20 years) and each year-month (for months in each year).

It generates a figure for each KO, containing subpanels of: A) heatmap of number of metagenomes for each month, B) 20 facets of AMG trend plots for each year-month, C) AMG trend plot for each month (aggregated across 20 years).

[script] 03.Reconstruct_vMAGs.step9.visualize_AMG_abundance.R

**10 Investigate AMG variation in each species and genus (non-singleton)**

The resulted file "Species_level_vOTU_AMG_variation.txt" contains the following information:

1) The number of non-singleton vOTU - excluding vOTUs with only one member

2) The number of non-singleton vOTU with any AMG KO hit(s)

3) The number of non-singleton vOTU with high AMG frequency (> 75%) in its all genomes - frequency here means the percentage of genomes containing KO among all the genomes

The resulted file "Species_level_vOTU_AMG_variation_statistics.txt" contains the following information:

1) The total number of vOTU and ko combinations (for vOTUs with any AMG KO hits(s))

2) The total number of vOTU and ko combinations (for vOTUs with any AMG KO hits(s)) of the 1st quartile size species-level vOTUs (The size range of 1st quartile size species-level vOTUs: 2-2)

3) The total number of vOTU and ko combinations (for vOTUs with any AMG KO hits(s)) of the 2nd quartile size species-level vOTUs (The size range of 2nd quartile size species-level vOTUs: 2-3)

4) The total number of vOTU and ko combinations (for vOTUs with any AMG KO hits(s)) of the 3rd quartile size species-level vOTUs (The size range of 3rd quartile size species-level vOTUs: 3-4)

5) The total number of vOTU and ko combinations (for vOTUs with any AMG KO hits(s)) of the 4th quartile size species-level vOTUs (The size range of 4th quartile size species-level vOTUs: 4-108)

The resulted file "Species_level_vOTU_AMG_variation_statistics.2.txt" [for KO distribution from the 1st quartile ko frequency (> 75%) and the 4th quartile in bin size] contains:

1) The total number of KO for AMG KO frequency > 75% and vOTU from the 4th quartile in bin size

2) The relative abundance of each KO (percentage of occurrence)

The resulted file "KO2ko_abun_n_mean_ko_freq.txt" [for KO distribution from the 1st quartile ko frequency (> 75%) and the 4th quartile in bin size] contains:

1) The relative abundance of each KO (percentage of occurrence)

2) The frequency of each KO (percentage of genomes containing this KO among all the genomes)

abundance and frequency were placed in two columns.

The resulted file "KO2dates_in_a_year.txt" [for KO distribution from the 1st quartile ko frequency (> 75%) and the 4th quartile in bin size] stores the presence/absence information of KO in a date of a year. 

[script] 03.Reconstruct_vMAGs.step10.investigate_AMG_variation_in_species.pl

