#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    import pandas as pd
    from pingouin import ttest
 
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)
    
# Aim: Conduct t test for pNpS results for 2000-2003 and 2016-2019

# Step 1 Conduct t test for pNpS results for four AMGs
pNpS_4AMGs = pd.read_csv('MetaPop.for_each_year/MetaPop/pNpS_result.for_2000-2003_and_2016-2019.for_four_AMGs.txt', sep = '\t', index_col = 0)
pNpS_4AMGs_dict = pNpS_4AMGs.to_dict() # col => row => value

viral_gene_4AMGs = list(pNpS_4AMGs.index.values)

year_2000_to_2003 = ['2000', '2001', '2002', '2003']
year_2016_to_2019 = ['2016', '2017', '2018', '2019']

genes_with_elevated_pNpS_4AMGs = {} # gene => [average_pNpS_2000_2003, average_pNpS_2016_2019, pval, cohend]
for gene in viral_gene_4AMGs:
    pNpS_2000_to_2003 = []
    pNpS_2016_to_2019 = []
    
    for year in year_2000_to_2003:
        if str(pNpS_4AMGs_dict[year][gene]) != 'NA' and str(pNpS_4AMGs_dict[year][gene]) != 'Inf' and str(pNpS_4AMGs_dict[year][gene]) != 'nan':
            pNpS_2000_to_2003.append(pNpS_4AMGs_dict[year][gene])
            
    for year in year_2016_to_2019:
        if str(pNpS_4AMGs_dict[year][gene]) != 'NA' and str(pNpS_4AMGs_dict[year][gene]) != 'Inf' and str(pNpS_4AMGs_dict[year][gene]) != 'nan':
            pNpS_2016_to_2019.append(pNpS_4AMGs_dict[year][gene])
            
    if len(pNpS_2000_to_2003) >= 3 and len(pNpS_2016_to_2019) >= 3: # At least containing 3 meaningful values        
        ttest_out = ttest(pNpS_2000_to_2003, pNpS_2016_to_2019)
        ttest_out = ttest_out.to_dict()
        pval = ttest_out['p-val']['T-test']
        cohend = ttest_out['cohen-d']['T-test']
        average_pNpS_2000_2003 = sum(pNpS_2000_to_2003) / len(pNpS_2000_to_2003)
        average_pNpS_2016_2019 = sum(pNpS_2016_to_2019) / len(pNpS_2016_to_2019)
        if pval != 'NaN' and cohend != 'NaN' and float(pval) < 0.05 and float(cohend) >= 0.7 and average_pNpS_2016_2019 > average_pNpS_2000_2003:
            genes_with_elevated_pNpS_4AMGs[gene] = [str(average_pNpS_2000_2003), str(average_pNpS_2016_2019), str(pval), str(cohend)]
            
f = open('MetaPop.for_each_year/MetaPop/genes_with_elevated_pNpS.for_2000-2003_and_2016-2019.for_four_AMGs.txt', 'w')
for gene in genes_with_elevated_pNpS_4AMGs:
    f.write(gene + '\t' + '\t'.join(genes_with_elevated_pNpS_4AMGs[gene]) + '\n')
f.close()   

# Step 2 Conduct t test for pNpS results for all viruses
pNpS = pd.read_csv('MetaPop.for_each_year/MetaPop/pNpS_result.for_2000-2003_and_2016-2019.txt', sep = '\t', index_col = 0)
pNpS_dict = pNpS.to_dict() # col => row => value

viral_gene = list(pNpS.index.values)

year_2000_to_2003 = ['2000', '2001', '2002', '2003']
year_2016_to_2019 = ['2016', '2017', '2018', '2019']

all_genes_with_elevated_pNpS = {} # gene => [average_pNpS_2000_2003, average_pNpS_2016_2019, pval, cohend]
for gene in viral_gene:
    pNpS_2000_to_2003 = []
    pNpS_2016_to_2019 = []
    
    for year in year_2000_to_2003:
        if str(pNpS_dict[year][gene]) != 'NA' and str(pNpS_dict[year][gene]) != 'Inf' and str(pNpS_dict[year][gene]) != 'nan':
            pNpS_2000_to_2003.append(pNpS_dict[year][gene])
            
    for year in year_2016_to_2019:
        if str(pNpS_dict[year][gene]) != 'NA' and str(pNpS_dict[year][gene]) != 'Inf' and str(pNpS_dict[year][gene]) != 'nan':
            pNpS_2016_to_2019.append(pNpS_dict[year][gene])
            
    if len(pNpS_2000_to_2003) >= 3 and len(pNpS_2016_to_2019) >= 3: # At least containing 3 meaningful values        
        ttest_out = ttest(pNpS_2000_to_2003, pNpS_2016_to_2019)
        ttest_out = ttest_out.to_dict()
        pval = ttest_out['p-val']['T-test']
        cohend = ttest_out['cohen-d']['T-test']
        average_pNpS_2000_2003 = sum(pNpS_2000_to_2003) / len(pNpS_2000_to_2003)
        average_pNpS_2016_2019 = sum(pNpS_2016_to_2019) / len(pNpS_2016_to_2019)
        if pval != 'NaN' and cohend != 'NaN' and float(pval) < 0.05 and float(cohend) >= 0.7 and average_pNpS_2016_2019 > average_pNpS_2000_2003:
            all_genes_with_elevated_pNpS[gene] = [str(average_pNpS_2000_2003), str(average_pNpS_2016_2019), str(pval), str(cohend)]
            
f = open('MetaPop.for_each_year/MetaPop/all_genes_with_elevated_pNpS.for_2000-2003_and_2016-2019.txt', 'w')
for gene in all_genes_with_elevated_pNpS:
    f.write(gene + '\t' + '\t'.join(all_genes_with_elevated_pNpS[gene]) + '\n')
f.close()  

