# Time-series analysis - Part 2 Each AMG analysis

(Get four important AMG-containing viral genome coverage statistics using time-series metagenomes across 20 years; the four important AMGs include: *psbA*, *pmoC*, *katG*, and *ahbD*)

**1 Get the psbA-containing viral genome coverage statistics**

We get the psbA-containing viral genome coverage statistics across all 465 metagenomes. The distribution of the viral genome should be >= 20 out of 465 metagenomes.

It will generate the following outputs:

1) PsbA_containing_viral_gn2distribution_n_cov_mean_n_tax_n_host_tax.txt

​    The distribution of psbA-containing viral genome number, mean coverage, viral taxonomy, and host taxonomy

2) PsbA_containing_viral_gn2month2cov.txt

​    The monthly coverage for psbA-containing viral genome

3) PsbA_containing_viral_gn2year_month2cov.txt

​    The year-month coverage for psbA-containing viral genome

4) PsbA_AMG_gene2month2cov_ratio.txt

​    The monthly coverage ratio for psbA AMG 

5) PsbA_AMG_gene2year_month2cov_ratio.txt

​    The year-month coverage ratio for psbA AMG 

[script] 08.Time_series_analysis.part2.each_AMG_analysis.step1.get_psbA-containing_viral_gn_coverage_statistics.pl

**2 Get the pmoC-containing viral genome coverage statistics**

The running process and outputs are similar to those described in Step 1 except that the distribution of the viral genome should be >= 5 out of 465 metagenomes.

[script] 08.Time_series_analysis.part2.each_AMG_analysis.step2.get_pmoC-containing_viral_gn_coverage_statistics.pl

**3 Get the katG-containing viral genome coverage statistics**

The running process and outputs are similar to those described in Step 1 except that the distribution of the viral genome should be >= 5 out of 465 metagenomes.

[script] 08.Time_series_analysis.part2.each_AMG_analysis.step3.get_katG-containing_viral_gn_coverage_statistics.pl

**4 Get the ahbD-containing viral genome coverage statistics**

The running process and outputs are similar to those described in Step 1 except that the distribution of the viral genome should be >= 20 out of 465 metagenomes.

[script] 08.Time_series_analysis.part2.each_AMG_analysis.step4.get_ahbD-containing_viral_gn_coverage_statistics.pl

