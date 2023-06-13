# TYMEFLIES_Viral
The repository stores scripts for the project of "TYMEFLIES_Viral" - Study of the viral population based on 20-year time series metagenome data from Lake Mendota, Madison, WI, US (metagenomes are obtained from lake water from pelagic integrated epilimnion zone).



Scripts (including some of the inputs/outputs) are placed in the following folders:

**1** [Process the datasets](https://github.com/AnantharamanLab/TYMEFLIES_Viral/blob/main/Processing_the_datasets): Copy fastq file, calculate fastq statistics, and get all metagenome assemblies cov state (or depth) files

**2** [Identify phages and find active prophages](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Identify_phages_and_find_active_prophages): Identify phages by VIBRANT, find active prophage by PropagAtE, and run CheckV to get phage scaffold quality

(This part is mainly based on the usage of software [VIBRANT](https://github.com/AnantharamanLab/VIBRANT), [PropagAtE](https://github.com/AnantharamanLab/PropagAtE), and [CheckV](https://bitbucket.org/berkeleylab/CheckV))

**3** [Reconstruct vMAGs](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Reconstruct_vMAGs): Recontruct vMAGs using vRhyme, get the best set of phage bins using stringent criteria, run CheckV to get phage vMAG quality, and summarize AMGs for all metagenomes

(This part is mainly based on the usage of software [vRhyme](https://github.com/AnantharamanLab/vRhyme); we also made a custom script to get the best set of phage bins)

**4** [Cluster phage genomes](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Cluster_phage_genomes): Cluster vMAGs into vOTUs at family-, genus-, and species-level

[The original phage vOTU clustering methods were adopted from two previously published papers: 1) [Nat Microbiol. 2021 Jul;6(7):960-970.](https://pubmed.ncbi.nlm.nih.gov/34168315/) 2) [Nucleic Acids Res. 2021 Jan 8;49(D1):D764-D775.](https://pubmed.ncbi.nlm.nih.gov/33137183/) While, since the large number of viral genomes (~1.3M genomes) in this study, we firstly clustered genomes into family- and genus-level vOTUs using [MCL-based method](https://github.com/snayfach/MGV/tree/master/aai_cluster) (we also modified the original python script within to reduce the RAM demand for our case), then used dRep to get species-level vOTUs within each genus.]

**5** [Taxonomic_classification](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Taxonomic_classification): Classify phage genomes using two methods: NCBI RefSeq viral protein searching and VOG HMM marker searching

(This part is mainly based on the method in [Nucleic Acids Res. 2021 Jan 8;49(D1):D764-D775.](https://pubmed.ncbi.nlm.nih.gov/33137183/))

**6** [Host prediction](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Host_prediction): Predict host using four approaches: 1) search sequence similarity to a microbial genome; 2) match to CRISPR spacers; 3) prophage scaffold search; 4) match to AMG (auxiliary metabolic gene)

(This part is mainly based on the method in [Nucleic Acids Res. 2021 Jan 8;49(D1):D764-D775.](https://pubmed.ncbi.nlm.nih.gov/33137183/))

**7** [Rscript for visualization](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Rscript_for_visualization): Rscripts for a variety of visualization works

**8** [Mapping metagenomic assemblies](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Mapping_metagenomic_assemblies): Map reads to metagenomic assemblies (the original scaffolds including both microbial and viral ones) to get MAG/virus abundance

**9** [Time series analysis - Part 1](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Time_series_analysis_part_1): Conduct AMG ratio and viral genome coverage analysis

**10** [Time series analysis - Part 2](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Time_series_analysis_part_2): Get four important AMG-containing viral genome coverage statistics

**11** [Time series analysis - Part 3](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Time_series_analysis_part_3): Get microdiversity analysis results



<br>

Database processing scripts are placed in the following folders:

**1** [Database IMGVR ](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Database_IMGVR): IMG/VR database v4.1 release Dec. 2022 (for [Cluster phage genomes](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Cluster_phage_genomes) and [Host prediction](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Host_prediction)) 

**2** [Database NCBI RefSeq viral](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Database_NCBI_RefSeq_viral):  NCBI RefSeq viral (2021-11-04 release) (for [Taxonomical classification](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Taxonomic_classification) and [Host prediction](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Host_prediction)) 

**3** [Database VOG97](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Database_VOG97):   VOG97 HMMs Release date Apr 19, 2021 (for [Taxonomical classification](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Taxonomic_classification)) 

**4** [Database TYMEFLIES MAGs](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Database_TYMEFLIES_MAGs): MAGs in IMG platform (for [Host prediction](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Host_prediction)) 

**5** [Database GEM](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Database_GEM): MAGs from GEM (A genomic catalog of Earth’s microbiomes) (for [Host prediction](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Host_prediction)) 

