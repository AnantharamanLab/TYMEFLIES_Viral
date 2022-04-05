# Time-series analysis

(Include AMG coverage ratio analysis and micro-/macrodiversity analysis for viruses by time-series metagenomes across 20 years)

**1** Get the information of high occurrence species (>= 40)

The occurrence info is based on the input file of "Species2occurrence_n_abundance.txt" (generated from XXX.pl in ). This occurrence is based on directly binned viral genomes in individual metagenomes. 

The information contains:1) Tax 2) Host tax 3) AMG KO within. This script also contains a filtering step: Only keep species that contain >= 10 genomes.

[script] 08.Time_series_analysis.step1.get_information_of_high_occurrence_species.pl