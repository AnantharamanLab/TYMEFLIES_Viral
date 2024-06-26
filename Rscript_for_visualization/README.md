# Rscripts for visualization

**1 Draw box plots for VIBRANT result**

The input file is "VIBRANT_result_summary.txt". The output figure is "VIBRANT_result_summary.png".

[script] 01.draw_boxplot_for_VIBRANT_result.R

**2 Draw box plots for phage length and completeness**

The input file is "phage_scaffold_and_vMAGs_length_and_completeness_info.txt" which was generated by "02.parse_to_get_phage_scaffold_and_vMAGs_length_and_completeness_info.pl". The output figure is "Phage_scaffold_or_bin_length_and_completeness.png".

[script] 02.parse_to_get_phage_scaffold_and_vMAGs_length_and_completeness_info.pl

[script] 03.draw_boxplot_for_phage_length_and_completeness.R

**3 Draw stack boxplot for scaffold binning state**

The input file is "Scaffold_to_binning_state_before_and_after_binning.txt" which was generated by "04.parse_to_get_scaffold_to_binning_state_before_and_after_binning.pl". The output figure is "Scaffold_binned_unbinned_percentage.png".

[script] 04.parse_to_get_scaffold_to_binning_state_before_and_after_binning.pl

[script] 05.draw_stack_boxplot_for_scaffold_binning_state.R

**4 Draw box plots for phage genome length for each CheckV quality category**

The input file is "Bin2checkv_quality_and_length.txt" which was generated by "06.parse_to_get_all_phage_geome_completeness_to_length_info.pl". The output figure is "Bin2checkv_quality_and_length.png".

[script] 06.parse_to_get_all_phage_geome_completeness_to_length_info.pl

[script] 07.draw_boxplot_for_phage_genome_length_for_each_checkv_quality_category.R

**5 Draw bar plots for bin member number (number of scaffolds in a bin) distribution**

The input file is "Bin_member_number_distribution.txt" which was generated by "08.parse_to_get_bin_member_number_distribution.pl". The output figure is "Bin_member_number_distribution.png".

Note: Only the bin member numbers with frequency > 1%  were shown in the resulted figure.

[script] 08.parse_to_get_bin_member_number_distribution.pl

[script] 09.draw_barplot_for_bin_member_number_distribution.R

**6 Draw bar plots for viral scaffold, vMAG, species, and genus number**

Summarize viral scaffold number, vMAG number, species-level vOTU number, and genus-level vOTU number. Note that viral genomes that are not from any genera were also added to the genus-level vOTU number. 

[script] 10.parse_to_get_scaffold_n_vMAG_n_species_n_genus_num.pl

[script] 11.draw_boxplot_for_Scaffold_vMAG_Species_Genus_num.R

**7 Draw species rarefaction curve**

The input file "Sample_num2species_num_for_rarefaction_curve.txt" was derived from the output of Step 7 of [Cluster phage genomes](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Cluster_phage_genomes) (nine replicates combined). It contains three columns: 1) Number of samples 2) Number of mean species found 3) Standard deviation of species found. 

The input file and output figure were also provided here:

Sample_num2species_num_for_rarefaction_curve.txt

Sample_num2species_num_for_rarefaction_curve.png

[script] 12.draw_rarefaction_curve.R

**8 Draw bar plots for seasonally viral family abundance, viral host family abundance, and MAG family abundance**

The input files and output figures were also provided here:

1) seasonally viral family abundance:

Season2family2abun.mdfed.txt

Season2family2abun.pdf

2) seasonally viral host family abundance

Season2host_family2abun.mdfed.txt

Season2host_family2abun.pdf

3) seasonally MAG family abundance

Family2season2abun_for_MAG.mdfed.txt

Season2family2abun_for_MAG.pdf

[script] 13.draw_bar_plots_for_seasonally_tax_n_host_tax_n_MAG_tax_distribution.R