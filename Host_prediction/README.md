# Host prediction

The method was based on *Nucleic Acids Res*. 2021 Jan 8;49(D1):D764-D775. doi: 10.1093/nar/gkaa946. 

(Link: https://academic.oup.com/nar/article/49/D1/D764/5952208).

It contains three approaches to find hosts: 

1) search sequence similarity to a microbial genome 

2) match to CRISPR spacers

3) prophage scaffold search



**1 Filter MAGs and find host by sequence similarity**

MAGs from TYMEFLIES and GEM databases are used as the reference to find virus-host association for the approach #1.  MAG contigs which were mainly viral (i.e. hit to a viral sequence at ≥90% identity over ≥ 50% of the host contig) were excluded as these can be incorrectly binned. Run blastn to find these potential viral contigs from all MAGs.

The queries are: TYMEFLIES and GEM MAGs (split into 20000-seqs fsa file individually for parallel running)

The blastn database (viral sequences) include: 

1) [IMG VR v3 all phages](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Database_IMGVR) 2) [NCBI RefSeq all viruses](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Database_NCBI_RefSeq_viral)
3) TYMEFLIES all phages 4) Lake Mendota time series (2008-2012) all phages (another dataset not shown)

After filtering viral contigs within MAG databases, virus-host predictions were then based on matches of ≥90% nucleotide identity covering ≥2 kb of both virus and (putative) host sequences. If there are multiple matches for a virus, 80% consensus at each taxonomic rank was used as the final host taxonomy.

[script] 06.Host_prediction.step1.run_blastn_to_filter_MAGs_and_identify_host_by_sequence_identity.pl

**2 Find CRISPR spacers from TYMEFLIES and GEM MAGs**

Two bioinformatics tools were used to identify CRISPRs from MAGs: 1) MinCED (https://github.com/ctSkennerton/minced) 2) PILER-CR (https://www.drive5.com/pilercr/) by using default settings. 

[script] 06.Host_prediction.step2.find_crispr_from_MAGs.pl

**3 Find MAG CRISPR spacer to virus matches by blastn**

Use blastn to search hits of spacer matches from all phage genome scaffolds for approach #2. Two categories of spacer matches were summarized:

1) have 0 or 1 mismatch over the entire spacer length (‘CRISPR (near)identical’) 
2) have ≥ 80% identity over the entire spacer length (‘CRISPR multiple partial’)

For each virus, a CRISPR-based prediction (at each taxonomic rank) was then derived based on an 80% consensus of all ‘CRISPR (near)identical’ predictions, or if no such prediction was available, an 80% consensus of ‘CRISPR multiple partial’ predictions. 

This is according to the methods in https://academic.oup.com/nar/article/49/D1/D764/5952208

[script] 06.Host_prediction.step3.run_blastn_to_find_matches.pl

**4 Find host for prophage**

For approach #3, since prophage scaffold were derived from its host, find host from TYMEFLIES MAGs that contained the scaffold where the prophage was located.

[script] 06.Host_prediction.step4.find_prophage_host.pl

