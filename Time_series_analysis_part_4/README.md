# Time-series analysis - Part 4 Virus and MAG taxa association analysis

(Get virus and MAG taxa association analysis results using time-series metagenomes across 20 years)

**1 Map all metagenomic reads to the collection of TYMEFLIES rep MAGs**

Use the metagenomic reads from 471 samples to map on to the collection of TYMEFLIES rep MAGs. The identity cutoff was preliminarily set to 90%.

[script] 08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step1.map_metagenomic_reads_to_the_collection_of_rep_MAG.pl

**2 get MAG abundance and family abundance**

Parse the mapping result from Step 1 to get MAG abundance and family abundance. The "--min-read-percent-identity" option of CoverM was set to 93% (due to that each MAG presents a species) to capture all the species level microdiversities. Process the MAG presence with breadth (>= 10%) cutoff.

[script] 08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step2.parse_to_get_MAG_abundance_and_family_abundance.py

**3** **Parse to get *psbA* virus and MAG abundance from 0 day of Early Summer for Cyanobiaceae**

(1) Parse to get the Cyanobiaceae MAG coverage from 0 day of Early Summer for each year

(2) Parse to get the *psbA*-containing Cyanobiaceae viral genome coverage from 0 day of Early Summer for each year

(3) Parse to get the no-*psbA*-containing Cyanobiaceae viral genome coverage from the 0 day of Early Summer for each year.

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step3.get_psbA_virus_and_MAG_abundance.for_Cyanobiaceae.py

**4** **Parse to get *psbA* virus and MAG abundance from 0 day of Early Summer for *Microcystis***

(1) Parse to get the *Microcystis* MAG coverage from 0 day of Early Summer for each year

(2) Parse to get the *psbA*-containing *Microcystis* viral genome coverage from 0 day of Early Summer for each year

(3) Parse to get the no-*psbA*-containing *Microcystis* viral genome coverage from the 0 day of Early Summer for each year.

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step4.get_psbA_virus_and_MAG_abundance.for_Microcystis.py

**5** **Parse to get *psbA* virus and MAG abundance from 0 day of Early Summer for *Planktothrix***

(1) Parse to get the *Planktothrix* MAG coverage from 0 day of Early Summer for each year

(2) Parse to get the *psbA*-containing *Planktothrix* viral genome coverage from 0 day of Early Summer for each year

(3) Parse to get the no-*psbA*-containing *Planktothrix* viral genome coverage from the 0 day of Early Summer for each year

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step5.get_psbA_virus_and_MAG_abundance.for_Planktothrix.py

**6** **Parse to get *psbA* virus and MAG abundance from 0 day of Early Summer for Other Cyanobacteria**

(1) Parse to get the Other Cyanobacteria MAG coverage from 0 day of Early Summer for each year

(2) Parse to get the *psbA*-containing Other Cyanobacteria viral genome coverage from 0 day of Early Summer for each year

(3) Parse to get the no-*psbA*-containing Other Cyanobacteria viral genome coverage from the 0 day of Early Summer for each year.

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step6.get_psbA_virus_and_MAG_abundance.for_Other_Cyanobacteria.py

**7** **Parse to get *pmoC* virus and MAG abundance from 0 day of Late Summer for *Methylocystis***

(1) Parse to get the *Methylocystis* MAG coverage from 0 day of Late Summer for each year

(2) Parse to get the *pmoC*-containing *Methylocystis* viral genome coverage from 0 day of Late Summer for each year

(3) Parse to get the no-*pmoC*-containing *Methylocystis* viral genome coverage from the 0 day of Late Summer for each year.

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step7.get_pmoC_virus_and_MAG_abundance.for_Methylocystis.py

**8** **Parse to get *pmoC* virus and MAG abundance from 0 day of Late Summer for UBA6136**

(1) Parse to get the UBA6136 MAG coverage from 0 day of Late Summer for each year

(2) Parse to get the *pmoC*-containing UBA6136 viral genome coverage from 0 day of Late Summer for each year

(3) Parse to get the no-*pmoC*-containing UBA6136 viral genome coverage from the 0 day of Late Summer for each year.

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step8.get_pmoC_virus_and_MAG_abundance.for_UBA6136.py

**9** **Parse to get *pmoC* virus and MAG abundance from 0 day of Late Summer for *Methylomonas***

(1) Parse to get the *Methylomonas* MAG coverage from 0 day of Late Summer for each year

(2) Parse to get the *pmoC*-containing *Methylomonas* viral genome coverage from 0 day of Late Summer for each year

(3) Parse to get the no-*pmoC*-containing *Methylomonas* viral genome coverage from the 0 day of Late Summer for each year.

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step9.get_pmoC_virus_and_MAG_abundance.for_Methylomonas.py

**10** **Parse to get *pmoC* virus and MAG abundance from 0 day of Late Summer for UBA10906**

(1) Parse to get the UBA10906 MAG coverage from 0 day of Late Summer for each year

(2) Parse to get the *pmoC*-containing UBA10906 viral genome coverage from 0 day of Late Summer for each year

(3) Parse to get the no-*pmoC*-containing UBA10906 viral genome coverage from the 0 day of Late Summer for each year.

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step10.get_pmoC_virus_and_MAG_abundance.for_UBA10906.py

**11** **Parse to get *katG* virus and MAG abundance from 0 day of Clearwater for *Planktophila***

(1) Parse to get the *Planktophila* MAG coverage from 0 day of Clearwater for each year

(2) Parse to get the *katG*-containing *Planktophila* viral genome coverage from 0 day of Clearwater for each year

(3) Parse to get the no-*katG*-containing *Planktophila* viral genome coverage from the 0 day of Clearwater for each year.

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step11.get_katG_virus_and_MAG_abundance.for_Planktophila.py

**12** **Parse to get *katG* virus and MAG abundance from 0 day of Clearwater for *Nanopelagicus***

(1) Parse to get the *Nanopelagicus* MAG coverage from 0 day of Clearwater for each year

(2) Parse to get the *katG*-containing *Nanopelagicus* viral genome coverage from 0 day of Clearwater for each year

(3) Parse to get the no-*katG*-containing *Nanopelagicus* viral genome coverage from the 0 day of Clearwater for each year.

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step12.get_katG_virus_and_MAG_abundance.for_Nanopelagicus.py

**13** **Parse to get *katG* virus and MAG abundance from 0 day of Clearwater for Other Nanopelagicales**

(1) Parse to get the Other Nanopelagicales MAG coverage from 0 day of Clearwater for each year

(2) Parse to get the *katG*-containing Other Nanopelagicales viral genome coverage from 0 day of Clearwater for each year

(3) Parse to get the no-*katG*-containing Other Nanopelagicales viral genome coverage from the 0 day of Clearwater for each year.

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step13.get_katG_virus_and_MAG_abundance.for_Other_Nanopelagicales.py

**14** **Parse to get *ahbD* virus abundance from 0 day of Clearwater**

Parse to get the *ahbD*-containing viral genome coverage from 0 day of Clearwater for each year

[script] 

08.Time_series_analysis.part4.virus_n_MAG_taxa_association.step14.get_ahbD_virus_abundance.py
