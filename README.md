# TYMEFLIES_Viral
The repository stores scripts for the project of "TYMEFLIES_Viral" - Study of the viral population based on 20-year time series metagenome data from Lake Mendota, Madison, Wisconsin, United States  (metagenomes are obtained from lake water from pelagic integrated epilimnion zone).

Scripts (including some of the inputs/outputs) are placed in the following folders:

1 [Process the datasets](https://github.com/AnantharamanLab/TYMEFLIES_Viral/blob/main/Processing_the_datasets): Copy fastq file, calculate fastq statistics, and get all metagenome assemblies cov state (or depth) files

2  [Identify phages and find active prophages](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Identify_phages_and_find_active_prophages): Identify phages by VIBRANT, find active prophage by PropagAtE, and run CheckV to get phage scaffold quality

3 [Reconstruct vMAGs](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Reconstruct_vMAGs): Recontruct vMAGs using vRhyme, get the best set of phage bins using stringent criteria, run CheckV to get phage vMAG quality, and summarize AMGs for all metagenomes

4 [Cluster vMAGs](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Cluster_vMAGs): Cluster vMAGs by dRep into vOTUs

5  [Taxonomic_classification](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Taxonomic_classification): Classify phage genomes using two methods: NCBI RefSeq viral protein searching and VOG HMM marker searching



Database processing scripts are placed in the following folders:

1 [Database IMGVR ](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Database_IMGVR): IMG_VR_2020-10-12_5.1 - IMG/VR v3 (for [Cluster vMAGs](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Cluster_vMAGs)) 

2 [Database NCBI RefSeq viral](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Database_NCBI_RefSeq_viral):  NCBI RefSeq viral (2021-11-04 release) (for [Taxonomical classification](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Taxonomic_classification)) 

3 [ Database VOG209](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Database_IMGVR):   VOG209 HMMs Release date Dec 07, 2021 (for [Taxonomical classification](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Taxonomic_classification)) 

