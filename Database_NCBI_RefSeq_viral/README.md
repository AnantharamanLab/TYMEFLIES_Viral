## NCBI RefSeq viral

1) NCBI RefSeq viral: Downloaded from https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/

The last modified date of this database is: 2021-11-04 (for [Taxonomical Classification](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Taxonomic_classification))

2) ICTV_Master_Species_List_2020.v1.txt is the 2021/5/10 release version (MSL #36)

The "viral.protein.ictv_8_rank_tax.txt" is the tax map for RefSeq viral proteins (not all proteins have a corresponding tax tag) in which the tax is in the style of ICTV (8 ranks in total)

3) "viral.protein.w_tax.faa" is a subset of "viral.protein.faa" with all proteins having a corresponding tax tag
viral.protein.w_tax.faa: 563742 sequences
viral.protein.faa 578023 sequences

4) The dmnd file (viral.protein.w_tax.dmnd) was made by diamond v0.9.14.115



[script] 01.grep_tax_info_from_gbff.pl

Parse tax information from gbff file (for all nucleotide sequences).

This tax is in the style of NCBI taxonomy.

[script] 02.grep_tax_info_from_gpff.pl

Parse tax information from gpff file (for all protein sequences).

This tax is in the style of NCBI taxonomy.

[script] 03.grep_NCBI_RefSeq_viral_proteins_w_tax.pl

Grep protein in "viral.protein.faa" to check whether a protein inside has its corresponding tax; if yes, we store them to "viral.protein.w_tax.faa".

[script] 04.reformat_NCBI_tax_to_ICTV_8_rank_tax.pl

Reformat NCBI tax of viruses to ICTV 8 rank tax (only contains "Realm Kingdom Phylum Class Order  Family  Genus Species")

[script] 05.make_diamond_blastp_db.sh

Make diamond db for "viral.protein.w_tax.faa".



[output] viral.protein.tax.txt.gz  - output of [script] 02.grep_tax_info_from_gpff.pl

[output] viral.protein.ictv_8_rank_tax.txt.gz  - output of [script] 04.reformat_NCBI_tax_to_ICTV_8_rank_tax.pl

