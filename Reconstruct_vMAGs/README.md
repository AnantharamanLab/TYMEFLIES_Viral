# Reconstruct vMAGs

**1 Run vRhyme to reconstruct vMAGs from phage scaffolds identified by VIBRANT**

Use default settings. vRhyme link: https://github.com/AnantharamanLab/vRhyme

[script] 09.run_vRhyme.pl

**2 Find the best bin set from vRhyme result**

Four criteria were used to refine the best bin set from vRhyme result:                     

1) Prophage identified by VIBRANT should be excluded from binning                        

2) Two or more lysogenic (non-prophage) phage scaffolds can not be in the same bin                          

3) Phage scaffolds identified by CheckV as "Complete" should be excluded from binning    

4) The maximum number of bin redundancy should be <= 1  

The resulted faa, ffn, fasta files for each vRhyme bin and vRhyme unbinned genome were provided

[script] 10.make_new_vRhyme_bins.pl

**3 Run CheckV for each phage genome**

Due to CheckV only allows for single-contig phage genome, we firstly link bin scaffolds with multiple Ns (1500 Ns as suggested by default) to make temporary single-contig phage genome. Then we run CheckV for temporary single-contig phage genome with default settings.

[script] 11.run_checkV_for_each_phage_genome.pl

[script] path/to/folder/vRhyme_v1.0.0/vRhyme/aux/link_bin_sequences.py (for linking bin scaffolds)

This script can be found in "https://github.com/AnantharamanLab/vRhyme/tree/master/vRhyme/aux"

**4 Summary CheckV results for each phage genome**

[script] 12.summarize_checkV_result_for_each_phage_genome.pl

**5 Summarize all AMGs from individual metagenomes** 

Step 1: Generate AMG table for all metagenomes. Use the AMG annotated by VIBRANT and change the protein according to the IMG ID and phage genome ID

Step 2: Add metagenome info and KEGG pathway info to the AMG table. The following information has been added:

Metagenome date and season ('date | date in the year | season')

Metabolisms (mutiple metabolisms separated by '|')

Pathways (mutiple pathways separated by '|')

[script] 13.summarize_AMG.step1.pl and 13.summarize_AMG.step2.pl



