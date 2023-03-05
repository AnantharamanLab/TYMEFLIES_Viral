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

**5 & 6 Summarize and filter all AMGs from individual metagenomes** 

Step 5: Generate AMG table for all metagenomes. Use the AMG annotated by VIBRANT and change the protein according to the IMG ID and phage genome ID

Step 6: Filter AMG and add metagenome info and KEGG pathway info to the AMG table. 

- Filtering step:

(1) Filter tail AMGs

Any one AMG that is placed in either ends of a scaffold or any n AMGs that are placed in either ends of a scaffold in tandem will be filtered

(2) Filter AMGs that have any v-scores (KEGG and Pfam v-scores) >= 1

In the original VIBRANT annotation result, we can parse the v-scores for individual AMG proteins. If any AMG has its any v-scores (KEGG, Pfam, and VOG v-scores) >= 1 (representing a viral-like nature), we will filter this AMG.

(3) Filter AMGs with flanking genes of v-scores < 0.02

For AMGs (or multiple AMGs placed in tandem) surrounded by all their flanking gene with v-scores < 0.02, it  indicates that these AMGs are surrounded by non-viral (cellular) genes. These AMGs will be filtered due to that they are likely to be non-viral (cellular) in origin too.

(4) Filter AMGs that have COG category as T or B

We used eggNOG-mapper v2  (in March 2023) to annotate all AMG proteins. We obtained the COG category information for all AMG proteins. If any AMG has its COG category assigned as T (Signal Transduction) or B (Chromatin Structure and dynamics), we  will filter this AMG.

[Result]

Total AMG number before filtering: 236,212

Total AMG number after filtering: 149,289

​          which means 63.2% retained

Number of AMGs to be filtered in each step:

*1* tail_AMG number is 46,438

*2* AMGs_w_high_v_scores number is 50,272 (using both KEGG v-score and Pfam v-score)

*3* AMGs_w_flanking_genes_of_low_v_score number is 96 (using KEGG v-score and considering two flanking genes at both sides)

*4* AMGs_belong_to_not_correct_COG number is 652

*Note that the above four collections of AMGs to be filtered have overlaps* 

- The following information has been added:

Metagenome date and season ('date | date in the year | season')

Metabolisms (mutiple metabolisms separated by '|')

Pathways (mutiple pathways separated by '|')

This script will use the file ["VIBRANT_KEGG_pathways_summary.tsv"](https://github.com/AnantharamanLab/VIBRANT/blob/master/files/VIBRANT_KEGG_pathways_summary.tsv) from VIBRANT.  

[script] 03.Reconstruct_vMAGs.step5.generate_AMG_summary_table.pl and 03.Reconstruct_vMAGs.step6.filter_AMG_and_add_metagenome_and_KEGG_info_to_AMG_summary_table.py

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

It generates a figure for each KO, containing subpanels of A) heatmap of the number of metagenomes for each month, B) 20 facets of AMG trend plots for each year-month, and C) AMG trend plot for each month (aggregated across 20 years).

[script] 03.Reconstruct_vMAGs.step9.visualize_AMG_abundance.R

**10 Investigate AMG variation in each species (non-singleton)**

The resulted file "Species_level_vOTU_AMG_variation.txt" contains the following information:

1) The number of non-singleton vOTU - excluding vOTUs with only one member (singleton vOTUs)

2) The number of non-singleton vOTU with any AMG KO hit(s)

3) The number of non-singleton vOTU with high AMG frequency (> 75%) in its all genomes - frequency here means the percentage of genomes containing the KO among all the genomes

The resulted file "Species_level_vOTU_AMG_variation_statistics.txt" contains the following information:

1) The total number of vOTU and KO combinations (for vOTUs with any AMG KO hits(s))

2) The total number of vOTU and KO combinations (for vOTUs with any AMG KO hits(s)) of the 1st quartile size species-level vOTUs (The size range of 1st quartile size species-level vOTUs: 2-2)

3) The total number of vOTU and KO combinations (for vOTUs with any AMG KO hits(s)) of the 2nd quartile size species-level vOTUs (The size range of 2nd quartile size species-level vOTUs: 2-3)

4) The total number of vOTU and KO combinations (for vOTUs with any AMG KO hits(s)) of the 3rd quartile size species-level vOTUs (The size range of 3rd quartile size species-level vOTUs: 3-4)

5) The total number of vOTU and KO combinations (for vOTUs with any AMG KO hits(s)) of the 4th quartile size species-level vOTUs (The size range of 4th quartile size species-level vOTUs: 4-108)

The resulted file "Species_level_vOTU_AMG_variation_statistics.2.txt" [for KO distribution from the 1st quartile KO frequency (> 75%) and the 4th quartile in bin size] contains:

1) The total number of KO from vOTU and KO combination collection of AMG KO frequency > 75% and vOTU from the 4th quartile in bin size

2) The relative abundance of each KO (relative abundance here means the percentage of the count of each KO divided by the total count of all KOs)

The resulted file "KO2ko_abun_n_mean_ko_freq.txt" [for KO distribution from the 1st quartile KO frequency (> 75%) and the 4th quartile in bin size] contains:

1) The relative abundance of each KO (relative abundance here means the percentage of the count of each KO divided by the total count of all KOs)

2) The mean frequency of each KO (frequency here means the percentage of genomes containing this KO among all the genomes)

Abundance and frequency were placed in two columns.

The resulted file "KO2dates_in_a_year.txt" [for KO distribution from the 1st quartile KO frequency (> 75%) and the 4th quartile in bin size] stores the presence/absence information of each KO on a date of a year. 

[script] 03.Reconstruct_vMAGs.step10.investigate_AMG_variation_in_species.pl

**11 Visualize AMG variation**

This script generates five plots:

1) General statistics of vOTU and KO combinations (bar plot)

[input] Species_level_vOTU_AMG_variation_statistics.table_1.txt

[output] Species_level_vOTU_AMG_variation_statistics.table_1.pdf

2) KO to vOTU tax abundance [bar plot, only nine KOs with the highest abundance (relative abundance of KO from vOTU and KO combination collection of the 1st quartile KO frequency (> 75%) and the 4th quartile in bin size) were chosen here]

[input] Species_level_vOTU_AMG_variation_statistics.table_2.txt

[output] Species_level_vOTU_AMG_variation_statistics.table_2.pdf

3) KO to vOTU host tax abundance [bar plot, only nine KOs with the highest abundance (relative abundance of KO from vOTU and KO combination collection of the 1st quartile KO frequency (> 75%) and the 4th quartile in bin size) were chosen here]

[input] Species_level_vOTU_AMG_variation_statistics.table_3.txt

[output] Species_level_vOTU_AMG_variation_statistics.table_3.pdf

4) AMG KO fraction (a.k.a, relative abundance) to the mean AMG KO frequency across all vOTUs (scatter plot)

[input] KO2ko_abun_n_mean_ko_freq.mdfed.txt

[output] KO2ko_abun_n_mean_ko_freq.pdf

5) Seasonal distribution of high occurrence KOs across metagenomes  

[input] KO2dates_in_a_year.mdfed.txt

[output] KO2dates_in_a_year.pdf

Only the "HighOccurrenceKO" (occurrence >= 400) and AMG KO fraction >= 1.25% ones were depicted in the figure.

The distribution of available metagenomes on the dates of a year was also depicted (the first row).

[script] 03.Reconstruct_vMAGs.step11.visualize_AMG_variation.R

All the inputs and outputs were also provided here.

**12 Find interested KO tax and host tax abundance info**

Find the KO to viral family abundance and KO to viral host family abundance.

Both KO abundance values were normalized by read number per metagenome, 100M reads per metagenome.

[script] 03.Reconstruct_vMAGs.step12.find_intereseted_KO_tax_n_host_tax_info.pl

**13 Visualize interested KO tax and host tax abundance**

This script generates two bar plots:

1) KO to viral family abundance for eight low diversity KOs

[input] KO2tax2abun_fraction.txt

[output] KO2tax2abun_fraction.8_low_diversity_KOs.pdf

2) KO to viral host family abundance for eight low diversity KOs

[input] KO2host_tax2abun_fraction.txt

[output] KO2host_tax2abun_fraction.8_low_diversity_KOs.pdf

[script] 03.Reconstruct_vMAGs.step13.visualize_interrested_KO_tax_n_host_tax_abundance.R

All the inputs and outputs were also provided here.

**14  Get KO to tax and host tax alpha diversity pattern**

The viral taxonomy for alpha diversity analysis was set to the family level. For each KO, we only randomly took 100 viral genomes with informative family-level taxonomy assigned.

The viral host taxonomy for alpha diversity analysis was set to the family level. For each KO, we only randomly took 50 viral genomes with informative family-level taxonomy assigned.

[script] 03.Reconstruct_vMAGs.step14.get_KO_to_tax_n_host_tax_alpha_diversity.pl

**15 Visualize KO to tax and host tax alpha diversity pattern**

This script generates four plots:

1) Scatter plot of KO alpha diversity (based on family) to KO occurrence using Shannon Index

[input] KO2family2viral_gn_num.txt and KO2occurrence_n_abundance.txt

[output] KO.family.shannon2occurrence.pdf

2) Scatter plot of KO alpha diversity (based on family) to KO occurrence using Simpson Index

[input] KO2family2viral_gn_num.txt and KO2occurrence_n_abundance.txt

[output] KO.family.simpson2occurrence.pdf

3) Scatter plot of KO alpha diversity (based on host family) to KO occurrence using Shannon Index

[input] KO2host_family2viral_gn_num.txt and KO2occurrence_n_abundance.txt

[output] KO.host_family.shannon2occurrence.pdf

4) Scatter plot of KO alpha diversity (based on host family) to KO occurrence using Simpson Index

[input] KO2host_family2viral_gn_num.txt and KO2occurrence_n_abundance.txt

[output] KO.host_family.simpson2occurrence.pdf

[script] 03.Reconstruct_vMAGs.step15.plot_KO_tax_n_host_tax_alpha_diversity.R

All the inputs and outputs were also provided here.

**16 Visualize AMG KO and host association**

This script depicts KO vs host abundance distribution across months, it generates four line charts:

1) KO vs host abundance distribution (K00507 and K01627 vs Cyanobacteria-Nostocaceae)

[input] Viral_host_association_table_1.txt

[output] Viral_host_association_table_1.pdf

2) KO vs host abundance distribution (K15895, K00798, and K01734 vs Bacteroidota-Flavobacteriaceae)

[input] Viral_host_association_table_2.txt

[output] Viral_host_association_table_2.pdf

3) KO vs host abundance distribution (K01666 vs Bacteroidota-UBA961)

[input] Viral_host_association_table_3.txt

[output] Viral_host_association_table_3.pdf

4) KO vs host abundance distribution (K02703 and K02706 vs Cyanobacteria-Cyanobiaceae)

[input] Viral_host_association_table_4.txt

[output] Viral_host_association_table_4.pdf

[script] 03.Reconstruct_vMAGs.step16.visualize_AMG_KO_host_association.R

All the inputs and outputs were also provided here.



