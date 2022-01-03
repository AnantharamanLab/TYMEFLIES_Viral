# VOG209

VOG209 HMMs
Release date    Dec 07, 2021
Release number  vog209
Data source     NCBI Refseq release 209

Number of HMMs in this release: 28406 (only the marker HMMs will be used)

(for [Taxonomical Classification](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Taxonomic_classification))

587 marker genes table ("VOG_marker_table.txt")
This table is derived from Supplementary Table S2 https://academic.oup.com/nar/article/49/D1/D764/5952208
Tax rank Allassoviricetes from "VOG_marker_table.txt" is not present in ICTV (MSL #36)
Tax rank Levivirales from "VOG_marker_table.txt" is not present in ICTV (MSL #36)
Tax rank Leviviridae from "VOG_marker_table.txt" is not present in ICTV (MSL #36)
Tax rank Luteoviridae from "VOG_marker_table.txt" is not present in ICTV (MSL #36)

This involved with 7 VOG markers which are not used in the downstream analysis:
VOG02065        sp|P03616|COAT_BPPRR Coat protein       Riboviria;Orthornavirae;Lenarviricota;Allassoviricetes;Levivirales;Leviviridae
VOG03747        sp|P07394|MATA_BPGA Maturation protein A        Riboviria;Orthornavirae;Lenarviricota;Allassoviricetes;Levivirales;Leviviridae
VOG03749        sp|P07393|RDRP_BPGA RNA-directed RNA polymerase beta chain      Riboviria;Orthornavirae;Lenarviricota;Allassoviricetes;Levivirales;Leviviridae
VOG00259        sp|P09504|P0_TYYVF Suppressor of silencing P0   Riboviria;Orthornavirae;Kitrinoviricota;Tolucaviricetes;Tolivirales;Luteoviridae
VOG00961        sp|P27578|CAPSD_BYDVN Major capsid protein      Riboviria;Orthornavirae;Kitrinoviricota;Tolucaviricetes;Tolivirales;Luteoviridae
VOG02650        sp|P09511|MVP_TYYVF Movement protein    Riboviria;Orthornavirae;Kitrinoviricota;Tolucaviricetes;Tolivirales;Luteoviridae
VOG05552        sp|P09514|MCAPS_TYYVF Minor capsid protein P3-RTD       Riboviria;Orthornavirae;Kitrinoviricota;Tolucaviricetes;Tolivirales;Luteoviridae

Riboviria;Orthornavirae;Lenarviricota;Allassoviricetes;Levivirales;Leviviridae was renamed to Riboviria;Orthornavirae;Lenarviricota;Leviviricetes;Norzivirales;Fiersviridae
according to ICTV2020, so the first three VOG tax were changed in the "VOG_marker_table.mdfed.txt"

Luteoviridae was abolished in 2020, so the last four VOG tax were renamed to Riboviria;Orthornavirae;NA;NA;NA;NA in the "VOG_marker_table.mdfed.txt"



[script] 01.check_if_all_VOG_marker_hmm_present.pl

Check if all VOG marker HMMs are present in release 209.

[script] 02.check_if_all_VOG_marker_hmm_tax_ranks_present_in_ICTV.pl

Check if all the VOG marker tax ranks are present in ICTV table, if not, we provide a modified table "VOG_marker_table.mdfed.txt" (7 places changed as indicated above).

[script] 03.hmmpress_in_batch.pl

Run hmmpress for all HMMs in batch.