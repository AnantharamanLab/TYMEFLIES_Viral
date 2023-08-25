# Metatranscriptome analysis

**1 Make the mapping reference for 24 hrs metatranscriptomic datasets (2015/8/20-21) based on corresponding metagenomes**

Use 3300042395 (2015/8/19) and 3300042470 (2015/8/22) as the mapping reference.

[script] 09.Metatranscriptome_analysis.step1.make_mapping_ref_based_on_corresponding_metaG.for_24hrs_metaT_dataset.py

**2 Make the mapping reference for 24 hrs metatranscriptomic datasets (2015/8/20-21) based on all virus species rep genomes** 

Use All_phage_species_rep_gn_containing_AMG.genes and AMG_counterpart_genes_and_flankings.fasta as the mapping reference; this mapping reference can also be used by the mapping procedure of 2020 metatranscriptomic datasets.

[script] 09.Metatranscriptome_analysis.step2.make_mapping_ref_based_on_all_phage_species_rep_gn.py

**3 Make the mapping reference for 2020 metatranscriptomic datasets based on corresponding metagenomes**

The metaT ID to corresponding metaG ID:

ME_2020_07_24_5m_C => 3300044631 (2018/7/24)

ME_2020_08_05_10M_B => 3300044608 (2018/8/6)

ME_2020_08_25_10M_B => 3300042355 (2018/8/25)

ME_2020_10_19_5M_B => 3300034116 (2018/10/24)



Use 3300044631 (2018/7/24) as the mapping reference for ME_2020_07_24_5m_C 

Use 3300044608 (2018/8/6) as the mapping reference for ME_2020_08_05_10M_B 

Use 3300042355 (2018/8/25) as the mapping reference for ME_2020_08_25_10M_B

Use 3300034116 (2018/10/24) as the mapping reference for ME_2020_10_19_5M_B

[script] 09.Metatranscriptome_analysis.step3.make_mapping_ref_based_on_corresponding_metaG.for_2020_metaT_dataset.py

**4 Map metaT to corresponding metaG for 2020 metaT datasets**

Use the above corresponding metaT to metaG datasets (as shown in Step 3) to do the mapping.

[script] 09.Metatranscriptome_analysis.step4.mapping_metaT_to_corresponding_metaG.for_2020_metaT_dataset.py

**5 Make the mapping idx for all_virus_species_rep_gn**

Use the mapping reference generated in Step 2 to make the mapping idx file.

[script] 09.Metatranscriptome_analysis.step5.create_mapping_idx_for_all_phage_species_rep_gn.py

**6 Map metaT to all_virus_species_rep_gn for 2020 metaT datasets**

[script] 09.Metatranscriptome_analysis.step6.mapping_metaT_to_all_phage_species_rep_gn.for_2020_metaT_dataset.py

**7 Get the viral gene expression results for 2020 metaT datasets using the corresponding metaG as the mapping reference**

This script parses the result of Step 4.

[script] 09.Metatranscriptome_analysis.step7.get_viral_gene_expression_results_based_on_corresponding_metaG.for_2020_metaT_dataset.py

**8 Get the viral gene expression results for 2020 metaT datasets using all virus species rep as the mapping reference**

This script parses the result of Step 6.

[script] 09.Metatranscriptome_analysis.step8.get_viral_gene_expression_results_based_on_all_phage_species_rep_gn.for_2020_metaT_dataset.py

**9 Map metaT to corresponding metaG for 24hrs metaT datasets (2015/8/20-21)**

The metaT ID to corresponding metaG ID:

diel_cycle_metaT => the_24hrs_metaT_ref

[script] 09.Metatranscriptome_analysis.step9.mapping_metaT_to_corresponding_metaG.for_24hrs_metaT_dataset.py

**10 Map metaT to all_virus_species_rep_gn for 24hrs metaT datasets (2015/8/20-21)**

[script]

09.Metatranscriptome_analysis.step10.mapping_metaT_to_all_phage_species_rep_gn.for_24hrs_metaT_dataset.py

**11 Get the viral gene expression results for 24hrs metaT datasets using the corresponding metaG as the mapping reference**

This script parses the result of Step 9.

[script] 

09.Metatranscriptome_analysis.step11.get_viral_gene_expression_results_based_on_corresponding_metaG.for_24hrs_metaT_dataset.py

**12 Get the viral gene expression results for 24 hrs metaT datasets using all virus species rep as the mapping reference**

This script parses the result of Step 10.

[script] 

09.Metatranscriptome_analysis.step12.get_viral_gene_expression_results_based_on_all_phage_species_rep_gn.for_24hrs_metaT_dataset.py











































