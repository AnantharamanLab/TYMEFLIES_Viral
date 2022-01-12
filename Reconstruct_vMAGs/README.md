# Reconstruct vMAGs

**1 Run vRhyme to reconstruct vMAGs from phage scaffolds identified by VIBRANT**

Use default settings. vRhyme link: https://github.com/AnantharamanLab/vRhyme

[script] 03.Reconstruct_vMAGs.step1.run_vRhyme.pl

**2 Find the best bin set from vRhyme result**

Four criteria were used to refine the best bin set from vRhyme result:                     

1) Prophage identified by VIBRANT should be excluded from binning                        

2) Two or more lysogenic (non-prophage) phage scaffolds can not be in the same bin                          

3) Phage scaffolds identified by CheckV as "Complete" should be excluded from binning    

4) The maximum number of bin redundancy should be <= 1  

The resulted faa, ffn, fasta files for each vRhyme bin and vRhyme unbinned genome were provided

[script] 03.Reconstruct_vMAGs.step2.make_new_vRhyme_bins.pl

**3 Run CheckV for each phage genome**

Due to CheckV only allows for single-contig phage genome, we firstly link bin scaffolds with multiple Ns (1500 Ns as suggested by default) to make temporary single-contig phage genome. Then we run CheckV for temporary single-contig phage genome with default settings.

[script] 03.Reconstruct_vMAGs.step3.run_checkV_for_each_phage_genome.pl

[script] path/to/folder/vRhyme_v1.0.0/vRhyme/aux/link_bin_sequences.py (for linking bin scaffolds)

This script can be found in "https://github.com/AnantharamanLab/vRhyme/tree/master/vRhyme/aux"

**4 Summary CheckV results for each phage genome**

[script] 13.Reconstruct_vMAGs.step4.summarize_checkV_result_for_each_phage_genome.pl

**5 Summarize all AMGs from individual metagenomes** 

Step 1: Generate AMG table for all metagenomes. Use the AMG annotated by VIBRANT and change the protein according to the IMG ID and phage genome ID

Step 2: Add metagenome info and KEGG pathway info to the AMG table. The following information has been added:

Metagenome date and season ('date | date in the year | season')

Metabolisms (mutiple metabolisms separated by '|')

Pathways (mutiple pathways separated by '|')

[script] 03.Reconstruct_vMAGs.step5.generate_AMG_summary_table.pl and 03.Reconstruct_vMAGs.step6.add_metagenome_and_KEGG_pathway_info_to_AMG_summary_table.pl

**6 Summarize AMG frequencies and abundance monthly and calculate AMG abundance for each year_month (for example, "2000-01")**

This script will generate following results:

1) KO2month2freq.txt

AMG frequencies (normalized by read number per metagenome and metagenome number per month)     for each month

2) KO2month2abun.txt

AMG abundances (normalized by read number per metagenome and metagenome number per month) for each month

3) KO2year_month2abun.txt

AMG abundances (normalized by read number per metagenome and metagenome number per year_month) for each year_month

4) KO2month_ko_details.txt

KO details: KO ID and KO detail information

5) Month2num_of_metagenomes.txt

the number of metagenomes in each month 

6) Year_month2num_of_metagenomes.txt

the number of metagenomes in each year_month 

[script] 03.Reconstruct_vMAGs.step7.summarize_AMG_frquency_and_abundance.pl

**7 Grep all AMG proteins and run against dbCAN2 and MEROPS databases to get annotated**

Get all AMG proteins to annotate by dbCAN2 and MEROPS databases.

This will add more columns to "AMG_summary.txt":

"CAZyme ID"

"Enzyme (CAZyme)"

"Substrate (CAZyme)"

"Class of Substrate (CAZyme)"

"MEROPS ID"

"MEROPS details"

[script] 03.Reconstruct_vMAGs.step8.annnoate_AMG_by_dbCAN2_and_MEROPS.pl



