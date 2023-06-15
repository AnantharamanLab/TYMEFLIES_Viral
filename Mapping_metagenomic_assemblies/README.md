# Mapping metagenomic assemblies

**1 Run Bowtie 2 to map reads onto metagenomic assemblies**

Map reads onto metagenomic assemblies (including both microbial and viral scaffolds, the original assemblies) for each metagenome

(Note that this step is computing-intensive)

[script] 07.Mapping_metagenomic_assemblies.step1.run_bowtie2_mapping_for_each_metagenome.pl

**2 Get MAG abundance pattern, including seasonly pattern (aggregated across 20 years), year-season pattern, and yearly pattern (for each year)**

These five files will be generated in the resulted folder "MAG_abundance":

Family2season2abun.txt - Family abundances (normalized by read number per metagenome, 100M reads per metagenome, and metagenome number per month) for each season

Family2year_season2abun.txt - Family abundances for each year-season (each season within a year, e.g., "2000-Spring"; normalizing steps were the same as the above)

Season2num_of_metagenomes.txt - number of metagenomes in each season

Year_season2num_of_metagenomes.txt - number of metagenomes in each year-season

Year2num_of_metagenomes.txt - number of metagenomes in each year

[script] 07.Mapping_metagenomic_assemblies.step2.get_MAG_abundance_pattern.pl

**3 Get the species-level vOTU occurrence and abundance distribution pattern**

The viral species (species-level vOTU) occurrence was based on counting the species member occurrence across all metagenomes.

The viral abundance was normalized by setting 100M reads per metagenome.

This script also contains a filtering step: Only keep species that contain >= 10 genomes.

[script] 07.Mapping_metagenomic_assemblies.step3.get_species_occurrence_n_abundance_distribution.pl

**4 Get the KO occurrence and abundance distribution pattern**

The KO abundance was the sum value of abundances of scaffolds containing this KO. KO abundance was also normalized by setting 100M reads per metagenome.

[script] 07.Mapping_metagenomic_assemblies.step4.get_KO_occurrence_n_abundance_distribution.pl

**5 Draw scatter plots for both species and KO occurrence and abundance by R**

Two PDFs were generated:

./AMG_analysis/Species2occurrence_n_abundance.pdf

./AMG_analysis/KO2occurrence_n_abundance.pdf

[script] 07.Mapping_metagenomic_assemblies.step5.draw_scatter_plots_for_species_n_KO_occurrence_n_abundance.R

**6 Get seasonal KO-carrying viral genome distribution**

Get the KO to season to KO-carrying genome distribution result 

KO-carrying genome distributing percentage = The number of IMGs that contain the ko / The number of IMGs that belong to the season × 100%

Note: Only high occurrence KO was include (occurrence >= 400)

[script]

07.Mapping_metagenomic_assemblies.step6.get_seasonal_KO_carrying_viral_genome_distribution.py

**7 Visualize the KO-carrying viral genome distribution in a heatmap**

A PDF figure was drawn:

./AMG_analysis/KO2season2KO_carrying_genome_distribution_percentage.pdf

[script] 07.Mapping_metagenomic_assemblies.step7.visualize_KO_carrying_viral_genome_distribution.R

**8 Get the seasonal viral tax (family level) and host tax (family level) distribution pattern**

These two files will be generated:

Season2family2abun.txt - The viral family abundance (normalized) for each season

Season2host_family2abun.txt - The viral host family abundance (it is the predicted host; still based on viral abundance, not directly the MAG abundance) for each season.

[script] 07.Mapping_metagenomic_assemblies.step6.get_seasonal_viral_tax_n_host_tax_distribution.pl