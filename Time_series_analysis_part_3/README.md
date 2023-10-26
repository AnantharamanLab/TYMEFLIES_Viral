# Time-series analysis - Part 3 Microdiversity analysis

(Get microdiversity analysis results using time-series metagenomes across 20 years)

**1 Get the AMG-containing viral genome microdiversity (nt diversity) statistics**

The distribution of viral genome should be >= 20 out of 471 metagenomes for calculating season and year_season distribution. 

One microdiversity parameter was obtained from this step: nucleotide (nt) diversity. 

[script] 08.Time_series_analysis.part3.microdiversity_analysis.step1.get_viral_gn_microdiversity_statistics.pl

**2 Get the AMG-containing viral genome SNP density statistics**

The distribution of viral genome should be >= 20 out of 471 metagenomes for calculating season and year_season distribution. 

One microdiversity parameter was obtained from this step: SNP density. 

[script] 08.Time_series_analysis.part3.microdiversity_analysis.step2.get_viral_gn_SNP_density_statistics.pl

**3 Conduct Spearman's correlation test between species cov with nt diversity and SNP density**

The script "calc_spearman_correlation.py" was used to conduct Spearman's correlation test between species cov with nt diversity and SNP density for both species containing four AMGs and all AMG-containing species.

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step3.conduct_correlation_analysis.pl

calc_spearman_correlation.py

**4 Get the viral species (four AMG containing viral species) pN/pS results**

Parse the "global_gene_microdiversity.tsv" file from MetaPop result folder (10.Microdiversity) to get pN/pS results for both all AMG-containing viral species and those four important AMG-containing viral species.

[script] 08.Time_series_analysis.part3.microdiversity_analysis.step4.get_viral_gn_pNpS_results.pl

**5 Conduct separate mapping for different seasons (summer and winter) for Fst calculations**

The script "filter_coverage_file.py" was used to filter sam files by setting a coverage cutoff.

The script "filter_bam_by_reference.py" was used to filter bam file by using the scaffold name of only viral species representatives.

Summer => Late Summer and Winter => Ice-on, separately.

[script] 08.Time_series_analysis.part3.microdiversity_analysis.step5.conduct_separate_mapping_for_summer_and_winter.pl

**6.1 Conduct Fst analysis for summer and winter sample groups by MetaPop**

Use the all summer and winter bam files as the input to conduct MetaPop analysis.

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step6.1.run_metapop.for_summer_vs_winter.pl

**6.2 Parse Fst result for winter group and summer group**

Viral scaffolds with Fst > 0.15 are identified as significantly differentiating between seasonal groups (also referred to as "fixed viral scaffolds").

The following requirements are used to obtain the selected genes located on the fixed Fst viral scaffolds:
(1) pi in winter group > pi in summer group (pi stands for nt diversity)
(2) gene N/S SNV ratio in winter group < gene N/S SNV ratio in summer group

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step6.2.parse_Fst_result_for_summer_and_winter.py

**7.1 Conduct Fst analysis for year 2000-2003 and year 2016-2019 by MetaPop**

Use the 2000-2003 and 2016-2019 merged bam files as the input to conduct MetaPop analysis.

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step7.1.run_metapop.for_2000-2003_vs_2016-2019.pl

**7.2 Parse Fst result for  year 2000-2003 and year 2016-2019**

Viral scaffolds with Fst > 0.15 are identified as significantly differentiating between seasonal groups (also referred to as "fixed viral scaffolds").

The following requirements are used to obtain the selected genes located on the fixed Fst viral scaffolds:
(1) pi in winter group > pi in summer group (pi stands for nt diversity)
(2) gene N/S SNV ratio in winter group < gene N/S SNV ratio in summer group

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step7.2.parse_Fst_result_for_2000-2003_and_2016-2019.py











**9 Get the viral species (four AMG containing viral species) pN/pS results for each year**

Parse the "global_gene_microdiversity.tsv" file from MetaPop result (for each year) folder (10.Microdiversity) to get pN/pS results for both all AMG-containing viral species and those four important AMG-containing viral species.

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step9.get_viral_gn_pNpS_results.for_each_year.pl

**10 Get the viral species (species rep gn) SNP allele frequencies for each year**

The reference allele was the majority allele in the year of 2018.

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step10.get_viral_gn_SNP_allele_freq.for_each_year.pl

**11 Parse the viral species (species rep gn) SNP allele frequencies for each year**

Get the linear regression patterns and Spearman's rank correlation patterns for each viral genome.

The script "calc_regression.py" was used to get the "reg_slope", "reg_yint", and "reg_rsquared" parameters.

The script "calc_spearman_correlation.py" was used to get the "spearman_corr" and "spearman_pval" parameters.

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step11.parse_viral_gn_SNP_allele_freq_pattern.for_each_year.pl

**12.1 Calculate gene coverage for AMG-containing viral gn (species representatives) for each year**

The following requirements were used:

(1) Only take sub-region of a given scaffold, cut the start and stop 150 bp regions

(2) If the gene length is less than 450 bp, we do not use it

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step12.1.calculate_gene_coverage.for_each_year.pl

**12.2 Calculate gene frequency for AMG-containing viral gn (species representatives) for each year**

The following requirements were used:

(1) Gene frequency was estimated as the coverage of each gene divided by the median coverage of all other genes in the genome

(2) The mean coverage of all other genes in the genome should be >= 5

(3) The gene number in a genome with a valid coverage (not "NA" and > 0) should be over 50% of the total gene number

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step12.2.calculate_gene_frequency.for_each_year.pl

**12.3 Get gene gain and loss for AMG-containing viral gn (species representatives) for each year**

Gene gain or loss was determined by comparing the mean gene frequency of year 2000-2003 and year 2016-2019 (the minimum difference of gene frequency to determine a gene gain or loss is 0.4).

[script] 08.Time_series_analysis.part3.microdiversity_analysis.step12.3.calculate_gene_gain_and_loss_result.for_each_year.pl

**13 Get the gene pN/pS results of 2000-2003 and 2016-2019 for viral species genomes**

Get the viral species genomes' gene pNpS results of 2000-2003 and 2016-2019 for comparing them to see if there are genes with elevated pNpS values

[script] 08.Time_series_analysis.part3.microdiversity_analysis.step13.compare_pNpS_results.for_2000-2003_and_2016-2019.pl

**14 Conduct t test for gene pN/pS results of 2000-2003 and 2016-2019 for viral species genomes**

Use the gene pN/pS results parsed from Step 13 to conduct t test to find genes with elevated pNpS values with significance support

[script] 08.Time_series_analysis.part3.microdiversity_analysis.step14.conduct_t_test_for_pNpS_results.for_2000-2003_and_2016-2019.py



















