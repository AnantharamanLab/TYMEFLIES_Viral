# TYMEFLIES_Viral
The repository stores scripts for the project of "TYMEFLIES_Viral" - Study of the viral population based on 20-year time series metagenome data from Lake Mendota, Madison, Wisconsin, United States  (metagenomes are obtained from lake water from pelagic integrated epilimnion zone)

Scripts (including some of the inputs/outputs) are placed in the following folders:

1 [Processing the datasets](https://github.com/AnantharamanLab/TYMEFLIES_Viral/blob/main/Processing_the_datasets): Copy fastq file, calculate fastq statistics, and get all metagenome assemblies cov state (or depth) files

2  [Identify phages and find active prophages](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Identify_phages_and_find_active_prophages): Identify phages by VIBRANT, find active prophage by PropagAtE, and run CheckV to get phage scaffold quality

3 [Reconstruct vMAGs](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Reconstruct_vMAGs): Recontruct vMAGs using vRhyme, get the best set of phage bins using stringent criteria, run CheckV to get phage vMAG quality, and summarize AMGs for all metagenomes

4 Cluster vMAGs: 

