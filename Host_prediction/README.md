# Host prediction

It contains three approaches to find hosts: 

1) iPHoP-based prediction

2) prophage scaffold search

3) match to AMG (auxiliary metabolic gene)



**1 Predict host using iPHoP**

The TYMEFLIES species representative MAGs (2855 MAGs total, dereplicated by dRep with 96% sequence identity cutoff) were added to the default iPHoP database “Sept_21_pub”.

[script] 06.Host_prediction.step1.run_iPHoP.py

**2 Parse iPHoP host prediction result**

Host-genome predictions were ultimately determined using the following guideline: the outcome (extracted from the result file: Host_prediction_to_genus_m90.csv) having the highest confidence score was designated as the definitive outcome in cases where multiple prediction outcomes for a single virus are available. 

[script] 06.Host_prediction.step2.parse_iPHoP_result.py

**3 Find host for prophage**

For approach #2, since prophage scaffold were derived from its host, find host from TYMEFLIES MAGs that contained the scaffold where the prophage was located.

[script] 06.Host_prediction.step3.find_prophage_host.pl

**4 Find host based on AMG**

For approach #3, We find viral host based on the AMG identity match between AMGs and microbial counterpart genes

[script] 06.Host_prediction.step4.find_host_based_on_AMG.pl

**5 Integrate all host prediction results**

Integrate all host prediction results (the host taxonomical results that are predicted from the above three methods should be at least down to the family level) and write down the final result. 



\[Method for assigning host taxonomy derived from vOTU host taxonomy\]: 

1) Get into each species cluster to see if any genomes have already got hits, then expand the host prediction to all the members within this species. 

2) Use prophage, AMG host prediction, and iPHoP results to get other vOTU (species) member host prediction.



The overlapped host taxonomies were solved down to the family level based on the following priority: 

1) provirus within a host genome; 2) AMG match to host genome; 3) iPHoP result; 4) derived from species host taxonomy

[script] 06.Host_prediction.step5.integrate_all_results.pl



