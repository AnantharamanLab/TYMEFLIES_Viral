#  **Identify phages and find active prophages** 

**1 Run VIBRANT to get viral (phage) scaffolds from metagenomes**

Use default settings, the minimum input scaffold length is 1000 bp. VIBRANT link: https://github.com/AnantharamanLab/VIBRANT

[script] 05.run_VIBRANT.pl

**2 Run PropagAtE to get active state of prophages**

Use default settings. PropagAtE link: https://github.com/AnantharamanLab/PropagAtE

[script] 06.run_PropagAtE.pl

**3 Run CheckV to get quality summary of all identified phage scaffolds**

Use default settings. CheckV link: https://bitbucket.org/berkeleylab/CheckV

[script] 07.run_checkV.pl

**4 Summarize VIBRANT result**

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









