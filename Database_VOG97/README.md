# VOG97

Release date: Apr 19, 2021

Release number: VOG97

Data source: NCBI Refseq release 97

Number of HMMs in this release: 25399 (only the marker HMMs will be used)

(for [Taxonomical classification](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Taxonomic_classification))

587 marker genes table ("VOG_marker_table.txt")
This table is derived from Supplementary Table S2 https://academic.oup.com/nar/article/49/D1/D764/5952208
Tax rank Allassoviricetes from "VOG_marker_table.txt" is not present in ICTV (MSL #37)
Tax rank Levivirales from "VOG_marker_table.txt" is not present in ICTV (MSL #37)
Tax rank Leviviridae from "VOG_marker_table.txt" is not present in ICTV (MSL #37)
Tax rank Luteoviridae from "VOG_marker_table.txt" is not present in ICTV (MSL #37)

This involved with 13 VOG markers which are not used in the downstream analysis:
VOG02065        sp|P03616|COAT_BPPRR Coat protein       Riboviria;Orthornavirae;Lenarviricota;Allassoviricetes;Levivirales;Leviviridae
VOG03747        sp|P07394|MATA_BPGA Maturation protein A        Riboviria;Orthornavirae;Lenarviricota;Allassoviricetes;Levivirales;Leviviridae
VOG03749        sp|P07393|RDRP_BPGA RNA-directed RNA polymerase beta chain      Riboviria;Orthornavirae;Lenarviricota;Allassoviricetes;Levivirales;Leviviridae
VOG00259        sp|P09504|P0_TYYVF Suppressor of silencing P0   Riboviria;Orthornavirae;Kitrinoviricota;Tolucaviricetes;Tolivirales;Luteoviridae
VOG00961        sp|P27578|CAPSD_BYDVN Major capsid protein      Riboviria;Orthornavirae;Kitrinoviricota;Tolucaviricetes;Tolivirales;Luteoviridae
VOG02650        sp|P09511|MVP_TYYVF Movement protein    Riboviria;Orthornavirae;Kitrinoviricota;Tolucaviricetes;Tolivirales;Luteoviridae
VOG05552        sp|P09514|MCAPS_TYYVF Minor capsid protein P3-RTD       Riboviria;Orthornavirae;Kitrinoviricota;Tolucaviricetes;Tolivirales;Luteoviridae
VOG00251        sp|P03724|GP14_BPT7 Internal virion protein gp14        Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Caudovirales;Autographiviridae
VOG01461        sp|Q05233|VG26_BPML5 Minor tail protein Gp26    Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Caudovirales;NA
VOG09528        sp|Q05219|VG13_BPML5 Gene 13 protein    Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Caudovirales;NA
VOG10818        sp|P20324|SCAF_BPT3 Capsid assembly scaffolding protein Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Caudovirales;Autographiviridae
VOG10819        sp|P03747|TUBE2_BPT7 Tail tubular protein gp12  Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Caudovirales;Autographiviridae
VOG11404        sp|P19195|Y9KD_BPBF2 Uncharacterized 9.2 kDa protein    Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;Caudovirales;Demerecviridae

Riboviria;Orthornavirae;Lenarviricota;Allassoviricetes;Levivirales;Leviviridae was renamed to Riboviria;Orthornavirae;Lenarviricota;Leviviricetes;Norzivirales;Fiersviridae
according to ICTV2021, so the first three VOG tax (1st-3rd) were changed in the "VOG_marker_table.mdfed.txt"

Luteoviridae was abolished in 2021, so the middle four VOG (4th-7th) tax were renamed to Riboviria;Orthornavirae;NA;NA;NA;NA in the "VOG_marker_table.mdfed.txt"

Caudovirales was abolished in 2021, so in the last six (8th-13th) tax, Caudovirales was renamed to "NA"

[script] 01.check_if_all_VOG_marker_hmm_present.pl

Check if all VOG marker HMMs are present in release 97.

[script] 02.check_if_all_VOG_marker_hmm_tax_ranks_present_in_ICTV.pl

Check if all the VOG marker tax ranks are present in ICTV table, if not, we provide a modified table "VOG_marker_table.mdfed.txt" (13 places changed as indicated above).

[script] 03.get_HMM_subset_and_hmmpress.pl

Get a subset of HMM of 587 VOG marker and run hmmpress.