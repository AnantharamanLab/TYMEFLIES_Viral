# Time-series analysis - Part 3 Microdiversity analysis

(Get microdiversity analysis results using time-series metagenomes across 20 years)

**1 Get the AMG-containing viral genome microdiversity (nt diversity) statistics**

The distribution of viral genome should be >= 20 out of 465 metagenomes for calculating month and year_month distribution. 

One microdiversity parameter was obtained from this step: nt diversity. 

[script] 08.Time_series_analysis.part3.microdiversity_analysis.step1.get_viral_gn_microdiversity_statistics.pl

**2 Get the AMG-containing viral genome SNP density statistics**

The distribution of viral genome should be >= 20 out of 465 metagenomes for calculating month and year_month distribution. 

One microdiversity parameter was obtained from this step: SNP density. 

[script] 08.Time_series_analysis.part3.microdiversity_analysis.step2.get_viral_gn_SNP_density_statistics.pl

**3 Conduct Spearman's correlation test between species cov with nt diversity and SNP density**

The script "calc_spearman_correlation.py" was used to conduct Spearman's correlation test.

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step3.conduct_correlation_analysis.pl

calc_spearman_correlation.py

**4 Get the viral species (four AMG containing viral species) pN/pS results**

Parse the "global_gene_microdiversity.tsv" file from MetaPop result folder (10.Microdiversity) to get pN/pS results for both all AMG-containing viral species and those four important AMG-containing viral species.

[script] 08.Time_series_analysis.part3.microdiversity_analysis.step4.get_viral_gn_pNpS_results.pl

**5 Conduct separate mapping for different seasons (summer and winter) for Fst calculations**

The script "filter_coverage_file.py" was used to filter sam files by setting a coverage cutoff.

The script "filter_bam_by_reference.py" was used to filter bam file by using the scaffold name of only viral species representatives

[script] 08.Time_series_analysis.part3.microdiversity_analysis.step5.conduct_separate_mapping_for_summer_and_winter.pl

**6.1 Conduct Fst analysis for summer and winter sample groups by inStrain_lite**

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step6.1.conduct_Fst_analysis_for_summer_and_winter_groups.pl

**6.2 Parse Fst result for winter group and summer group**

The following requirements are used to parse the result from Step 6.1:
(1) Fst >= 0.5
(2) pi in winter group > pi in summer group (pi stands for nt diversity)
(3) gene N/S SNV ratio in winter group < gene N/S SNV ratio in summer group
(4) gene coverages in winter group and summer group both > 5×

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step6.2.parse_Fst_result_for_summer_and_winter_groups.pl

**7.1 Conduct Fst analysis for year 2000-2003 and year 2016-2019 by inStrain_lite**

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step7.1.conduct_Fst_analysis_for_2000-2003_and_2016-2019.pl

**7.2 Parse Fst result for  year 2000-2003 and year 2016-2019**

The following requirements are used to parse the result from Step 7.1:
(1) Fst >= 0.5
(2) pi in year 2000-2003 > pi in year 2016-2019
(3) gene N/S SNV ratio in year 2000-2003 < gene N/S SNV ratio in year 2016-2019
(4) gene coverages in year 2000-2003 and 2016-2019 both > 5×

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step7.2.parse_Fst_result_for_2000-2003_and_2016-2019.pl

**8 Parse the viral species Fst from "MetaPop.for_each_year" folder**

Two outputs were generated:

(1) Viral_scf2year_pair_between_2000_2003_and_2016_2019_period2fst_high.txt

The high Fst (>=0.5) result between 2000-2003 and 2016-2019 period

(2) Viral_scf2year_pair_between_2000_2003_and_2016_2019_period2fst_high.for_four_AMGs.txt

The high Fst  (>=0.5) result between 2000-2003 and 2016-2019 period for viral species containing four AMGs

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step8.parse_viral_gn_Fst_by_MetaPop.pl

**9 Get the viral species (four AMG containing viral species) pN/pS results for each year**

Parse the "global_gene_microdiversity.tsv" file from MetaPop result (for each year) folder (10.Microdiversity) to get pN/pS results for both all AMG-containing viral species and those four important AMG-containing viral species.

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step9.get_viral_gn_pNpS_results.for_each_year.pl

**10 Get the viral species (species rep gn) SNP allele frequencies for each year**

The reference allele was the majority allele in the year of 2018.

[script] 

08.Time_series_analysis.part3.microdiversity_analysis.step10.get_viral_gn_SNP_allele_freq.for_each_year.pl

**11 Parse the viral species (species rep gn) SNP allele frequencies for each year**

Get the linear regression patterns for each viral genome.

The script "calc_regression.py" was used to get the "reg_slope", "reg_yint", and "reg_rsquared" parameters.

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

08.Time_series_analysis.part3.microdiversity_analysis.step12.3.calculate_gene_gain_and_loss_result.for_each_year.pl













