# Host prediction

It contains three approaches to find hosts: 

1) iPHoP-based prediction

2) prophage scaffold search

4) match to AMG (auxiliary metabolic gene)



**1 Predict host using iPHoP**

The TYMEFLIES species representative MAGs (2855 MAGs total, dereplicated by dRep54 with 96% sequence identity cutoff) were added to the default iPHoP database “Sept_21_pub”.

[script] 06.Host_prediction.step1.run_iPHoP.py

**2 Parse iPHoP host prediction result**

The host prediction to genome results (based on host-based tools, including “blast”, “CRISPR”, and “iPHoP-RF” results) were finally assigned with the following rules: 1) If  “blast” or “CRISPR” results were obtained for one virus genome, the result with the highest confidence score was assigned as the final result. 2) If only “iPHoP-RF” results were obtained for one virus genome, the result with the highest confidence score was assigned as the final result. Note: we did not include the phage-based method (RaFAH) result. 

[script] 06.Host_prediction.step2.parse_iPHoP_result.py

**3 Find host for prophage**

For approach #2, since prophage scaffold were derived from its host, find host from TYMEFLIES MAGs that contained the scaffold where the prophage was located.

[script] 06.Host_prediction.step3.find_prophage_host.pl

**4 Find host based on AMG**

For approach #3, We find viral host based on the AMG identity match between AMGs and microbial counterpart genes

[script] 06.Host_prediction.step4.find_host_based_on_AMG.pl

**5 Integrate all host prediction results**

Integrate all host prediction results and write down the final result. 



\[Method for assigning host taxonomy derived from vOTU host taxonomy\]: 

1) Get into each species cluster to see if any genomes have already got hits, then expand the host prediction to all the members within this species cluster. 

2) Only use prophage, AMG host prediction, and "blast" or "CRISPR"-based iPHoP results to get other vOTU member host prediction.



The overlapped host taxonomies were solved based on the following priority: 

1) provirus within a host genome; 2) “blast”-based iPHoP result; 3) “CRISPR”-based iPHoP result; 4) AMG match to host genome; 5) “iPHoP-RF” result; 6) derived from species host taxonomy  

This priority order was mainly based on the description in the reference (*Nucleic Acids Res*. 2021 Jan 8;49(D1):D764-D775. doi: 10.1093/nar/gkaa946).

[script] 06.Host_prediction.step5.integrate_all_results.pl



