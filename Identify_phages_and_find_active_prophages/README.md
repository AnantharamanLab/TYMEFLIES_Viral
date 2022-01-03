#  **Identify phages and find active prophages** 

**1 Run VIBRANT to get viral (phage) scaffolds from metagenomes**

Use default settings, the minimum input scaffold length is 1000 bp. VIBRANT link: https://github.com/AnantharamanLab/VIBRANT

[script] 05.run_VIBRANT.pl

**2 Parse VIBRANT result folder to do modifications**

1) Prophage results from VIBRANT need to be refined:

Cutoff for two prophages in the same saffold to be combined and this scaffold is a entire phage scaffold: 
1. protein distance <= 20 proteins or nt distance <= 20 kb                                                
2. a fragment gap in between both prophages (for example, prophages are fragment_1 and fragment_3)        
3. host region (front, end, and gap regions counted together) <= 30%                                      

Cutoff for three or more prophages:                                                                       
1. protein distance <= 20 proteins or nt distance <= 20 kb                                                
2. host region (front, end, and gap regions counted together) <= 30%    

2) Set the minimal phage scaffold length to be >= 2000 bp. These phages and relevant information will be written in a separated folder named "VIBRANT_$img_id.a.v2.min$min_length". 

[script] Side04.amend_script_to_VIBRANT.pl

**3 Run PropagAtE to get active state of prophages**

Use default settings. PropagAtE link: https://github.com/AnantharamanLab/PropagAtE

[script] 06.run_PropagAtE.pl

**4 Run CheckV to get quality summary of all identified phage scaffolds**

Use default settings. CheckV link: https://bitbucket.org/berkeleylab/CheckV

[script] 07.run_checkV.pl

**5 Summarize VIBRANT result**

Summarize VIBRANT result for each IMG metagenome, here only count phage scaffold >= 2000 bp. 

Counted items:

Total scaffold num

Scaffold num over min2000

Phage num in total

Lytic phage num

Lysogenic phage (excluding prophage) num

Prophage num

Complete phage num

[script] 08.summarize_VIBRANT_result.v2.pl









