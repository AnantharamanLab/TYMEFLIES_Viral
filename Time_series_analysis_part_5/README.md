# Time-series analysis - Part 5 Virus, host, and environmental parameter association analysis

(Conduct the virus, host, and environmental parameter association analysis using time-series metagenomes across 20 years)

**1 Get the year seasonal abundance results of psbA AMG, cyanobacteria virus, and cyanobacteria**

Get the year-season abundance for:
	(1) psbA AMG cov (parse intermediate results from "AMG_gene2IMG2cov_ratio.txt" and "Viral_gn2IMG2cov_norm_filtered.txt")
	(2) psbA-containing Cyanobacteria virus abundance
	(3) no-psbA-containing Cyanobacteria virus abundance  
	(4) Cyanobacteria abundance

[script] 08.Time_series_analysis.part5.virus_n_host_n_env_association.step1.get_year_seasonal_abundance.py

**2 Calculate the duration days with cyanobacteria and cyanobacteria virus abundance maintained at > 20% of the peak abundance**

Conduct correlation analysis between each pair of virus/host and environmental parameters. We use the new method - try to find whether there are correlations between the number of days maintained at > 20% of the peak abundance (both cyanobacteria and cyanobacteria virus) and the mean environmental parameters during the summer seasons (both Early Summer and Late Summer).  Duration days were determined by computing the abundance profile using an interpolation function for each year, employing intervals of five days. The interpolation function was applied to abundance data within the time range of -45 to 160 days (since the start of Early Summer for each year). The years that could not meet the full time range were excluded from the correlation analysis. 

[script] 08.Time_series_analysis.part5.virus_n_host_n_env_association.step2.conduct_correlation_analysis.calculate_duration_date_for_Cyanobacteria.py

**3 Conduct Spearman's ranking correlation test for the duration days and the mean environmental parameters during the summer seasons**

The mean value of each parameter for every season was calculated by averaging all the measurements within that specific season (Some missing values were denoted as “NA”).  

The environmental parameters were obtained from both Early Summer and Late Summer, and their significance was weighted based on the dates of each summer season. The Spearman’s rank correlation test with Fisher z-transformation was employed, and the *p*-value was provided. The input table for biological entity vs environmental parameter correlation analysis was provided in this folder as "virus_n_host_n_env_association_input_table.for_Cyanobacteria.txt". The correlation result was visualized by R script (R script and the input file provided in this folder) using heatmap. 

[script] 08.Time_series_analysis.part5.virus_n_host_n_env_association.step3.conduct_correlation_analysis.conduct_spearman_ranking_correlation_test_for_Cyanobacteria.py

01.draw_heatmap_for_virus_n_host_n_env_association.R

