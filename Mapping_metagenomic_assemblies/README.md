# Mapping metagenomic assemblies

**1 Run Bowtie 2 to map reads onto metagenomic assemblies**

Map reads onto metagenomic assemblies (including both microbial and viral scaffolds, the original assemblies) for each metagenome

(Note that this step is computing-intensive)

[script] 07.Mapping_metagenomic_assemblies.step1.run_bowtie2_mapping_for_each_metagenome.pl

**2 Get MAG abundance pattern, both monthly pattern (aggregated across 20 years) and yearly pattern (for months in each year)**

These four files will be generated in the resulted folder "MAG_abundance":

Family2month2abun.txt - Family abundances (normalized by read number per metagenome, 100M reads per metagenome, and metagenome number per month) for each month

Family2year_month2abun.txt - Family abundances for each year-month (each month within a year, e.g., "2000-09"; normalizing steps were the same as the above)

Month2num_of_metagenomes.txt - number of metagenomes in each month

Year_month2num_of_metagenomes.txt - number of metagenomes in each year-month

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

[script] 07.Mapping_metagenomic_assemblies.step5.draw_scatter_plots_for_species_n_KO_occurrence_n_abundance.R

**6 Get the monthly viral tax (family level) and host tax (family level) distribution pattern**

These two files will be generated:

Month2family2abun.txt - The viral family abundance (normalized) for each month

Month2host_family2abun.txt - The viral host family abundance (it is the predicted host; still based on viral abundance, not directly the MAG abundance) for each month.

[script] 07.Mapping_metagenomic_assemblies.step6.get_monthly_viral_tax_n_host_tax_distribution.pl