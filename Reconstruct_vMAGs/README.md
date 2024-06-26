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

Any AMG placed at either ends of a scaffold or any number of AMGs placed at either ends of a scaffold in tandem were filtered.   

(2) Filter AMGs that have any v-scores (KEGG and Pfam v-scores) >= 1

The KEGG v-scores and Pfam v-scores for individual AMGs were parsed out from the VIBRANT result. If any AMG had its any v-scores (KEGG or Pfam v-scores) ≥ 1 (representing a viral-like nature), this AMG was filtered.   

(3) Filter AMGs with flanking genes of v-scores (only KEGG v-scores) < 0.25

For any AMG (or multiple AMGs placed in tandem) surrounded by all their flanking genes with v-scores < 0.25 (only KEGG v-scores considered here), it indicated that these AMGs were surrounded by non-viral (cellular) genes. These AMGs were filtered due to that they were likely to be non-viral (cellular) in origin  

(4) Filter AMGs that have COG category as T or B

The eggNOG-mapper v2 (in March 2023) was used to annotate all AMG proteins to get COG category assignment. If any AMG had its COG category assigned as “T” (Signal Transduction) or “B” (Chromatin Structure and dynamics), this AMG was then filtered.

[Result]

Total AMG number before filtering: 238,050

Total AMG number after filtering: 150,458

​          which means 63.2% retained

Number of AMGs to be filtered in each step:

*1* tail_AMG number is 46,765

*2* AMGs_w_high_v_scores number is 50,712 (using both KEGG v-score and Pfam v-score)

*3* AMGs_w_flanking_genes_of_low_v_score number is 97 (using KEGG v-score and considering two flanking genes at both sides)

*4* AMGs_belong_to_not_correct_COG number is 657

*Note that the above four collections of AMGs to be filtered have overlaps* 

- The following information has been added:

Metagenome date and season ('date | date in the year | season')

Metabolisms (multiple metabolisms separated by '|')

Pathways (multiple pathways separated by '|')

This script will use the file ["VIBRANT_KEGG_pathways_summary.tsv"](https://github.com/AnantharamanLab/VIBRANT/blob/master/files/VIBRANT_KEGG_pathways_summary.tsv) from VIBRANT.  

[script] 03.Reconstruct_vMAGs.step5.generate_AMG_summary_table.pl and 03.Reconstruct_vMAGs.step6.filter_AMG_and_add_metagenome_and_KEGG_info_to_AMG_summary_table.py

Note: The AMG summary table were further manually filtered (removing *queC/E/F* from the final AMG list due to that they are for the biosynthesis of queuosine which is against host restriction system). The final AMG number is 143,751.

**7 Get all virus normalized abundance**

First, store all the normalized depths for individual scaffolds from the mapping result of all-scaffolds-mapping. Then, get the mean scaffold normalized depth for all scaffolds within a viral genome.

[script] 03.Reconstruct_vMAGs.step7.get_all_virus_abundance.py

**8 Summarize AMG frequencies and abundance monthly (seasonal) and calculate AMG abundance for each year-month (year-season) (for example, "2000-01")**

This script will generate following results:

1) KO2month2freq.txt

AMG frequencies (normalized by read number per metagenome, 100M reads per metagenome, and metagenome number per month) for each month

2) KO2month2abun.txt

AMG abundances (normalized by read number per metagenome, 100M reads per metagenome, and metagenome number per month) for each month

3) KO2year_month2abun.txt

AMG abundances (normalized by read number per metagenome, 100M reads per metagenome, and metagenome number per year-month) for each year-month

4) KO2month_ko_details.txt

Two columns containing KO ID and KO detail for all KOs that appear in all metagenomes

5) Month2num_of_metagenomes.txt

the number of metagenomes in each month 

6) Year_month2num_of_metagenomes.txt

the number of metagenomes in each year_month 

[script] 03.Reconstruct_vMAGs.step8.summarize_AMG_frequency_and_abundance.pl

03.Reconstruct_vMAGs.step8.summarize_AMG_frequency_and_abundance.v2.pl (This script will generate similar results but with "seasonal" or "year-season" numbers)

**9 Grep all AMG proteins and run against dbCAN2 and MEROPS databases to get annotated**

Get all AMG proteins to annotate by dbCAN2 and MEROPS databases.

This will add more columns to "AMG_summary.txt":

"CAZyme ID"

"Enzyme (CAZyme)"

"Substrate (CAZyme)"

"Class of Substrate (CAZyme)"

"MEROPS ID"

"MEROPS details"

[script] 03.Reconstruct_vMAGs.step9.annnoate_AMG_by_dbCAN2_and_MEROPS.pl

**10 Make AMG trend R plots for each month (season) and each year-month (year-season)**

Use the results generated from Step 8 ("KO2month2abun.txt", "KO2year_month2abun.txt", and "KO2month_ko_details.txt") to draw AMG trend R plots (AMG KO abundance) for each month (aggregated across 20 years) and each year-month (for months in each year).

It generates a figure for each KO, containing subpanels of A) heatmap of the number of metagenomes for each month, B) 20 facets of AMG trend plots for each year-month, and C) AMG trend plot for each month (aggregated across 20 years).

[script] 03.Reconstruct_vMAGs.step10.visualize_AMG_abundance.R

03.Reconstruct_vMAGs.step10.visualize_AMG_abundance.v2.R (This script will generate similar results but with "seasonal" or "year-season" figures)

**11 Investigate AMG variation in each species (non-singleton)**

The resulted file "Species_level_vOTU_AMG_variation.txt" contains the following information:

1) The number of non-singleton vOTU - excluding vOTUs with only one member (singleton vOTUs)

2) The number of non-singleton vOTU with any AMG KO hit(s)

3) The number of non-singleton vOTU with high AMG frequency (> 75%) in its all genomes - frequency here means the percentage of genomes containing the KO among all the genomes

The resulted file "Species_level_vOTU_AMG_variation_statistics.txt" contains the following information:

1) The total number of vOTU and KO combinations (for vOTUs with any AMG KO hits(s))

2) The total number of vOTU and KO combinations (for vOTUs with any AMG KO hits(s)) of the 1st quartile size species-level vOTUs (The size range of 1st quartile size species-level vOTUs: 2-2)

3) The total number of vOTU and KO combinations (for vOTUs with any AMG KO hits(s)) of the 2nd quartile size species-level vOTUs (The size range of 2nd quartile size species-level vOTUs: 2-3)

4) The total number of vOTU and KO combinations (for vOTUs with any AMG KO hits(s)) of the 3rd quartile size species-level vOTUs (The size range of 3rd quartile size species-level vOTUs: 3-4)

5) The total number of vOTU and KO combinations (for vOTUs with any AMG KO hits(s)) of the 4th quartile size species-level vOTUs (The size range of 4th quartile size species-level vOTUs: 4-110)

The resulted file "Species_level_vOTU_AMG_variation_statistics.2.txt" [for KO distribution from the 1st quartile KO frequency (> 75%) and the 4th quartile in bin size] contains:

1) The total number of KO from vOTU and KO combination collection of AMG KO frequency > 75% and vOTU from the 4th quartile in bin size

2) The relative abundance of each KO (relative abundance here means the percentage of the count of each KO divided by the total count of all KOs)

The resulted file "KO2ko_abun_n_mean_ko_freq.txt" [for KO distribution from the 1st quartile KO frequency (> 75%) and the 4th quartile in bin size] contains:

1) The relative abundance of each KO (relative abundance here means the percentage of the count of each KO divided by the total count of all KOs)

2) The mean frequency of each KO (frequency here means the percentage of genomes containing this KO among all the genomes)

Abundance and frequency were placed in two columns.

[script] 03.Reconstruct_vMAGs.step11.investigate_AMG_variation_in_species.pl

**12 Visualize AMG variation**

This script generates two plots:

1) General statistics of vOTU and KO combinations (bar plot)

[input] Species_level_vOTU_AMG_variation_statistics.table_1.txt

[output] Species_level_vOTU_AMG_variation_statistics.table_1.pdf

2) AMG KO fraction (a.k.a, relative abundance) to the mean AMG KO frequency across all vOTUs (scatter plot)

[input] KO2ko_abun_n_mean_ko_freq.mdfed.txt

[output] KO2ko_abun_n_mean_ko_freq.pdf

[script] 03.Reconstruct_vMAGs.step12.visualize_AMG_variation.R

All the inputs and outputs were also provided here.

**13 Find interested KO tax and host tax abundance info**

Find the KO to viral family abundance and KO to viral host family abundance.

Both KO abundance values were normalized by read number per metagenome, 100M reads per metagenome.

[script] 03.Reconstruct_vMAGs.step13.find_intereseted_KO_tax_n_host_tax_info.pl

**14  Get KO to tax and host tax alpha diversity pattern**

The viral taxonomy for alpha diversity analysis was set to the family level. For each KO, we only randomly took 100 viral genomes with informative family-level taxonomy assigned.

The viral host taxonomy for alpha diversity analysis was set to the family level. For each KO, we only randomly took 25 viral genomes with informative family-level taxonomy assigned.

[script] 03.Reconstruct_vMAGs.step14.get_KO_to_tax_n_host_tax_alpha_diversity.pl

**15 Visualize interested KO tax and host tax abundance**

This script generates two bar plots:

1) KO to viral family abundance for eight low diversity KOs

[input] KO2tax2abun_fraction.txt

[output] KO2tax2abun_fraction.8_low_diversity_KOs.pdf

2) KO to viral host family abundance for eight low diversity KOs

[input] KO2host_tax2abun_fraction.txt

[output] KO2host_tax2abun_fraction.8_low_diversity_KOs.pdf

[script] 03.Reconstruct_vMAGs.step15.visualize_interrested_KO_tax_n_host_tax_abundance.R

All the inputs and outputs were also provided here.

**16 Visualize KO to tax and host tax alpha diversity pattern**

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

[script] 03.Reconstruct_vMAGs.step16.plot_KO_tax_n_host_tax_alpha_diversity.R

All the inputs and outputs were also provided here.

