# Cluster phage genomes

**1 Cluster phage genomes by shared protein an AAI**

Firstly, perform all-vs-all DIAMOND BLASTP for all phage genome proteins. Then, calculate the shared proteins  and the average AAI of shared proteins between genome pairs.

This script is modified from "amino_acid_identity.py" (https://github.com/snayfach/MGV/tree/master/aai_cluster). It was custom to use low RAM.

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