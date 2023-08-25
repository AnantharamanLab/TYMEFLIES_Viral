# Time-series analysis - Part 4 Virus and MAG taxa association analysis

(Get microdiversity analysis results using time-series metagenomes across 20 years)

**1 Map all metagenomic reads to the collection of TYMEFLIES rep MAGs**

Use the metagenomic reads from 471 samples to map on to the collection of TYMEFLIES rep MAGs. The identity cutoff was preliminarily set to 90%.

[script] 08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step1.map_metagenomic_reads_to_the_collection_of_rep_MAG.pl

**2 get MAG abundance and family abundance**

Parse the mapping result from Step 1 to get MAG abundance and family abundance. The "--min-read-percent-identity" option of CoverM was set to 93% (due to that each MAG presents a species) to capture all the species level microdiversities. Process the MAG presence with both coverage (>= 0.33) and breadth (>= 50%) cutoffs and also filter the scaffold abundance by removing the highest 5% and lowest 5% fractions.

[script] 08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step2.parse_to_get_MAG_abundance_and_family_abundance.py

**3** **Parse to get *psbA* virus and MAG abundance from 0 day of Early Summer**

(1) Parse to get the Cyanobiaceae MAG coverage from 0 day of Early Summer for each year

(2) Parse to get the *psbA*-containing viral genome coverage from 0 day of Early Summer for each year

(3) Parse to get the non-*psbA*-containing viral genome with Cyanobiaceae host coverage from the 0 day of Early Summer for each year.

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step3.get_psbA_virus_and_MAG_abundance.py

**4 Parse to get *pmoC* virus and MAG abundance from 0 day of Late Summer**

(1) Parse to get the Methanotroph MAG coverage from 0 day of Late Summer for each year

(2) Parse to get the *pmoC*-containing viral genome coverage from 0 day of Late Summer for each year

(3) Parse to get the non-*pmoC*-containing viral genome with Methanotroph host coverage from the 0 day of Late Summer for each year.

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step4.get_pmoC_virus_and_MAG_abundance.py



**The descriptions for 08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step5.get_katG_virus_abundance.py and** 

**08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step6.get_ahbD_virus_abundance.py are left. Finish them in the future!**

**8/24**





















