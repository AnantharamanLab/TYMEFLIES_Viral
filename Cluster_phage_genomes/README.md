# Cluster phage genomes

**1 Cluster viral genomes by shared protein an AAI**

Firstly, perform all-vs-all DIAMOND BLASTP for all viral genome proteins. Then, calculate the shared proteins and the average AAI of shared proteins between genome pairs.

This script is modified from "amino_acid_identity.py" (https://github.com/snayfach/MGV/tree/master/aai_cluster). It was custom for low RAM.

[script] 04.Cluster_phage_genomes.step1.cluster_phage_genomes_by_shared_protein_and_AAI.pl

**2 Parse shared protein and AAI result to make MCL input**

Use "filter_aai.py" (https://github.com/snayfach/MGV/tree/master/aai_cluster) to make MCL inputs (filter edges) for both genus clusters and family clusters.

[script] 04.Cluster_phage_genomes.step2.parse_to_get_MCL_input.pl

**3 Perform MCL analysis**

The settings are according to Section "Perform MCL-based clustering" (https://github.com/snayfach/MGV/tree/master/aai_cluster) for both genus clustering and family clustering.

[script] 04.Cluster_phage_genomes.step3.perform_MCL.pl

**4 Compare the genus cluster and family cluster result**

Get the following results:

Total number of viral genomes (in genus clusters)

Total number of viral genomes (in family clusters)

Total number of genus cluster

Total number of family cluster

Total number of singleton genus cluster (only one genome in a genus cluster)

The average percentage of genomes from each genus cluster fall into the same family cluster

The average percentage of genomes from each genus cluster fall into the same family cluster (excluding singleton genus clusters)

[script] 04.Cluster_phage_genomes.step4.compare_genus_cluster_and_family_cluster.pl

**5 Get vOTU (species level) from each genus cluster**

Use dRep to cluster viral genomes for individual genus clusters.

The final species-level vOTUs contain:

1) species-level vOTUs from each genus cluster

2) singleton genus genome

3) viral genome that is not included in any genera

[script] 04.Cluster_phage_genomes.step5.get_vOTU_from_each_genus_cluster.pl

**6 Make the input file for genus rarefaction curve**

We firstly picked a random sample from all 471 samples, and counted the number of genera that were found in this sample. Then, added samples gradually to all the 471 samples, made a table storing the sample number to genus number. We replicated this 49 times, and generated 49 resulting files.

[script] 04.Cluster_phage_genomes.step6.make_rarefaction_curve_for_genus.pl

**7 Make the input file for species rarefaction curve**

We firstly picked a random sample from all 471 samples, and counted the number of species that were found in this sample. Then, added samples gradually to all the 471 samples, made a table storing the sample number to species number. We replicated this 9 times, and generated 9 resulting files.

Visualization of species rarefaction curve can be found in Step 12 of [Rscript for visualization](https://github.com/AnantharamanLab/TYMEFLIES_Viral/tree/main/Rscript_for_visualization).

[script] 04.Cluster_phage_genomes.step7.make_rarefaction_curve_for_species.pl



