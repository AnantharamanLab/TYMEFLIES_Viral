#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    from pathlib import Path 
    from collections import defaultdict
    warnings.filterwarnings('ignore')
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
# Aim: Get the viral gene expression results for 2020 metaT datasets
# The metaT ID to corresponding metaG ID:
# ME_2020_07_24_5m_C => 3300044631 (2018/7/24)
# ME_2020_08_05_10M_B => 3300044608 (2018/8/6)
# ME_2020_08_25_10M_B => 3300042355 (2018/8/25)
# ME_2020_10_19_5M_B => 3300034116 (2018/10/24)

# Step 1 Store the gene function map (only contain genes that are from viral genomes)
Gene_fun_map = {} # gene => [viral_gene, annotation]
with open('/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_2020_metaT_ref_gene_function_map.based_on_MetaG.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        tmp = line.split('\t')
        if tmp[1] != 'NA':
            gene, viral_gene, annotation = tmp[0], tmp[1], tmp[2]
            Gene_fun_map[gene] = [viral_gene, annotation]
lines.close()            
            
# Step 2 Parse the coverage result in 'ME_2020_07_24_5m_C___3300044631_mapping_result_dir' folder   
Coverage_result_ME_2020_07_24_5m_C___3300044631 = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_07_24_5m_C___3300044631_mapping_result_dir/all_coverm_raw_result.txt'
Gene_expression_ME_2020_07_24_5m_C___3300044631_dict = {} # gene => [gene, viral_gene, annotation, gene_expression]
with open(Coverage_result_ME_2020_07_24_5m_C___3300044631, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Ga'):
            tmp = line.split('\t')
            gene, gene_expression = tmp[0], tmp[1]
            if gene in Gene_fun_map: # If the gene is from any viral genomes
                viral_gene, annotation = Gene_fun_map[gene][0], Gene_fun_map[gene][1]
                Gene_expression_ME_2020_07_24_5m_C___3300044631_dict[gene] = [gene, viral_gene, annotation, gene_expression]
lines.close()

## Step 2.1 Write down the gene expression dict
header = 'gene\tviral_gene\tannotation\tgene_expression(rpkm)\n'
f = open('/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_07_24_5m_C___3300044631_mapping_result_dir/viral_gene_expression_result.txt','w')
f.write(header)
for gene in Gene_expression_ME_2020_07_24_5m_C___3300044631_dict:
    line = '\t'.join(Gene_expression_ME_2020_07_24_5m_C___3300044631_dict[gene]) + '\n'
    f.write(line)
f.close()   

## Step 2.2 Get the gene expression result for viral genomes that have at least 50% genes with positive gene expression values
Gene_expression_ME_2020_07_24_5m_C___3300044631_viral_genomes_selected2gene_expression_list = defaultdict(lambda: [[], []]) # viral_gn => [[AMGs],[non-AMGs]]
for gene in Gene_expression_ME_2020_07_24_5m_C___3300044631_dict:   
    viral_gene, annotation, gene_expression = Gene_expression_ME_2020_07_24_5m_C___3300044631_dict[gene][1], Gene_expression_ME_2020_07_24_5m_C___3300044631_dict[gene][2], Gene_expression_ME_2020_07_24_5m_C___3300044631_dict[gene][3]
    viral_gn = viral_gene.split('__')[0] + '__' + viral_gene.split('__')[1]
    if annotation == 'NA':
        Gene_expression_ME_2020_07_24_5m_C___3300044631_viral_genomes_selected2gene_expression_list[viral_gn][1].append(gene_expression)
    elif annotation != 'NA':
        Gene_expression_ME_2020_07_24_5m_C___3300044631_viral_genomes_selected2gene_expression_list[viral_gn][0].append(gene_expression)
        
viral_genomes_not_selected = [] # The viral genomes not pass with positive expression
for viral_gn in Gene_expression_ME_2020_07_24_5m_C___3300044631_viral_genomes_selected2gene_expression_list:
    total_gene_no = 0
    total_gene_w_positive_expression_no = 0
    total_gene_no = len(Gene_expression_ME_2020_07_24_5m_C___3300044631_viral_genomes_selected2gene_expression_list[viral_gn][0]) + len(Gene_expression_ME_2020_07_24_5m_C___3300044631_viral_genomes_selected2gene_expression_list[viral_gn][1])
    for gene_expression in Gene_expression_ME_2020_07_24_5m_C___3300044631_viral_genomes_selected2gene_expression_list[viral_gn][0]:
        if float(gene_expression) > 0:
            total_gene_w_positive_expression_no = total_gene_w_positive_expression_no + 1
    for gene_expression in Gene_expression_ME_2020_07_24_5m_C___3300044631_viral_genomes_selected2gene_expression_list[viral_gn][1]:
        if float(gene_expression) > 0:
            total_gene_w_positive_expression_no = total_gene_w_positive_expression_no + 1
    if float(total_gene_w_positive_expression_no / total_gene_no) < 0.5:
        viral_genomes_not_selected.append(viral_gn)
        
for viral_gn in viral_genomes_not_selected:
    if viral_gn in Gene_expression_ME_2020_07_24_5m_C___3300044631_viral_genomes_selected2gene_expression_list:
        del Gene_expression_ME_2020_07_24_5m_C___3300044631_viral_genomes_selected2gene_expression_list[viral_gn]
        
## Step 2.3 Get the number of viral genomes (with positive expression) and the number of viral genomes that have positive AMG gene expression
no_of_viral_gns_w_positive_expression = 0
no_of_viral_gns_having_positive_AMG_gene_expression = 0
no_of_viral_gns_w_positive_expression = len(Gene_expression_ME_2020_07_24_5m_C___3300044631_viral_genomes_selected2gene_expression_list)
for viral_gn in Gene_expression_ME_2020_07_24_5m_C___3300044631_viral_genomes_selected2gene_expression_list:
    AMG_gene_expression_list = Gene_expression_ME_2020_07_24_5m_C___3300044631_viral_genomes_selected2gene_expression_list[viral_gn][0]
    if any(AMG_gene_expression_list): # The any function returns True if any element of an iterable is True, and False otherwise.
        no_of_viral_gns_having_positive_AMG_gene_expression = no_of_viral_gns_having_positive_AMG_gene_expression + 1

## Step 2.4 Write the result of Step 2.3        
f = open('/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_07_24_5m_C___3300044631_mapping_result_dir/no_of_viral_gn_w_positive_expression_and_AMG_gene_expression.txt','w')
header2 = 'No. of viral genomes with positive expression\tNo. of viral genomes having positive AMG gene expression\n'
f.write(header2)
f.write(f'{no_of_viral_gns_w_positive_expression}\t{no_of_viral_gns_having_positive_AMG_gene_expression}\n')
f.close()


# Step 3 Parse the coverage result in 'ME_2020_08_05_10M_B___3300044608_mapping_result_dir' folder   
Coverage_result_ME_2020_08_05_10M_B___3300044608 = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_08_05_10M_B___3300044608_mapping_result_dir/all_coverm_raw_result.txt'
Gene_expression_ME_2020_08_05_10M_B___3300044608_dict = {} # gene => [gene, viral_gene, annotation, gene_expression]
with open(Coverage_result_ME_2020_08_05_10M_B___3300044608, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Ga'):
            tmp = line.split('\t')
            gene, gene_expression = tmp[0], tmp[1]
            if gene in Gene_fun_map: # If the gene is from any viral genomes
                viral_gene, annotation = Gene_fun_map[gene][0], Gene_fun_map[gene][1]
                Gene_expression_ME_2020_08_05_10M_B___3300044608_dict[gene] = [gene, viral_gene, annotation, gene_expression]
lines.close()

## Step 3.1 Write down the gene expression dict
header3 = 'gene\tviral_gene\tannotation\tgene_expression(rpkm)\n'
f = open('/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_08_05_10M_B___3300044608_mapping_result_dir/viral_gene_expression_result.txt','w')
f.write(header3)
for gene in Gene_expression_ME_2020_08_05_10M_B___3300044608_dict:
    line = '\t'.join(Gene_expression_ME_2020_08_05_10M_B___3300044608_dict[gene]) + '\n'
    f.write(line)
f.close()   

## Step 3.2 Get the gene expression result for viral genomes that have at least 50% genes with positive gene expression values
Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list = defaultdict(lambda: [[], []]) # viral_gn => [[AMGs],[non-AMGs]]
for gene in Gene_expression_ME_2020_08_05_10M_B___3300044608_dict:   
    viral_gene, annotation, gene_expression = Gene_expression_ME_2020_08_05_10M_B___3300044608_dict[gene][1], Gene_expression_ME_2020_08_05_10M_B___3300044608_dict[gene][2], Gene_expression_ME_2020_08_05_10M_B___3300044608_dict[gene][3]
    viral_gn = viral_gene.split('__')[0] + '__' + viral_gene.split('__')[1]
    if annotation == 'NA':
        Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list[viral_gn][1].append(gene_expression)
    elif annotation != 'NA':
        Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list[viral_gn][0].append(gene_expression)
        
viral_genomes_not_selected = set() # The viral genomes not pass with positive expression and do not contain any AMGs
for viral_gn in Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list:
    total_gene_no = 0
    total_gene_w_positive_expression_no = 0
    total_gene_no = len(Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list[viral_gn][0]) + len(Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list[viral_gn][1])
    for gene_expression in Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list[viral_gn][0]:
        if float(gene_expression) > 0:
            total_gene_w_positive_expression_no = total_gene_w_positive_expression_no + 1
    for gene_expression in Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list[viral_gn][1]:
        if float(gene_expression) > 0:
            total_gene_w_positive_expression_no = total_gene_w_positive_expression_no + 1
    if float(total_gene_w_positive_expression_no / total_gene_no) < 0.5:
        viral_genomes_not_selected.add(viral_gn)
    if len(Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list[viral_gn][0]) == 0: # Store the viral genomes do not contain any AMGs
        viral_genomes_not_selected.add(viral_gn)
        
for viral_gn in viral_genomes_not_selected:
    if viral_gn in Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list:
        del Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list[viral_gn]
        
## Step 3.3 Get the number of viral genomes (with positive expression) and the number of viral genomes that have positive AMG gene expression
no_of_viral_gns_w_positive_expression = 0
no_of_viral_gns_having_positive_AMG_gene_expression = 0
no_of_viral_gns_w_positive_expression = len(Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list)
for viral_gn in Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list:
    AMG_gene_expression_list = Gene_expression_ME_2020_08_05_10M_B___3300044608_viral_genomes_selected2gene_expression_list[viral_gn][0]
    if any(AMG_gene_expression_list): # The any function returns True if any element of an iterable is True, and False otherwise.
        no_of_viral_gns_having_positive_AMG_gene_expression = no_of_viral_gns_having_positive_AMG_gene_expression + 1

## Step 3.4 Write the result of Step 3.3        
f = open('/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_08_05_10M_B___3300044608_mapping_result_dir/no_of_viral_gn_w_positive_expression_and_AMG_gene_expression.txt','w')
header4 = 'No. of viral genomes with positive expression\tNo. of viral genomes having positive AMG gene expression\n'
f.write(header4)
f.write(f'{no_of_viral_gns_w_positive_expression}\t{no_of_viral_gns_having_positive_AMG_gene_expression}\n')
f.close()


# Step 4 Parse the coverage result in 'ME_2020_08_25_10M_B___3300042355_mapping_result_dir' folder   
Coverage_result_ME_2020_08_25_10M_B___3300042355 = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_08_25_10M_B___3300042355_mapping_result_dir/all_coverm_raw_result.txt'
Gene_expression_ME_2020_08_25_10M_B___3300042355_dict = {} # gene => [gene, viral_gene, annotation, gene_expression]
with open(Coverage_result_ME_2020_08_25_10M_B___3300042355, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Ga'):
            tmp = line.split('\t')
            gene, gene_expression = tmp[0], tmp[1]
            if gene in Gene_fun_map: # If the gene is from any viral genomes
                viral_gene, annotation = Gene_fun_map[gene][0], Gene_fun_map[gene][1]
                Gene_expression_ME_2020_08_25_10M_B___3300042355_dict[gene] = [gene, viral_gene, annotation, gene_expression]
lines.close()

## Step 4.1 Write down the gene expression dict
header5 = 'gene\tviral_gene\tannotation\tgene_expression(rpkm)\n'
f = open('/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_08_25_10M_B___3300042355_mapping_result_dir/viral_gene_expression_result.txt','w')
f.write(header5)
for gene in Gene_expression_ME_2020_08_25_10M_B___3300042355_dict:
    line = '\t'.join(Gene_expression_ME_2020_08_25_10M_B___3300042355_dict[gene]) + '\n'
    f.write(line)
f.close()   

## Step 4.2 Get the gene expression result for viral genomes that have at least 50% genes with positive gene expression values
Gene_expression_ME_2020_08_25_10M_B___3300042355_viral_genomes_selected2gene_expression_list = defaultdict(lambda: [[], []]) # viral_gn => [[AMGs],[non-AMGs]]
for gene in Gene_expression_ME_2020_08_25_10M_B___3300042355_dict:   
    viral_gene, annotation, gene_expression = Gene_expression_ME_2020_08_25_10M_B___3300042355_dict[gene][1], Gene_expression_ME_2020_08_25_10M_B___3300042355_dict[gene][2], Gene_expression_ME_2020_08_25_10M_B___3300042355_dict[gene][3]
    viral_gn = viral_gene.split('__')[0] + '__' + viral_gene.split('__')[1]
    if annotation == 'NA':
        Gene_expression_ME_2020_08_25_10M_B___3300042355_viral_genomes_selected2gene_expression_list[viral_gn][1].append(gene_expression)
    elif annotation != 'NA':
        Gene_expression_ME_2020_08_25_10M_B___3300042355_viral_genomes_selected2gene_expression_list[viral_gn][0].append(gene_expression)
        
viral_genomes_not_selected = [] # The viral genomes not pass with positive expression
for viral_gn in Gene_expression_ME_2020_08_25_10M_B___3300042355_viral_genomes_selected2gene_expression_list:
    total_gene_no = 0
    total_gene_w_positive_expression_no = 0
    total_gene_no = len(Gene_expression_ME_2020_08_25_10M_B___3300042355_viral_genomes_selected2gene_expression_list[viral_gn][0]) + len(Gene_expression_ME_2020_08_25_10M_B___3300042355_viral_genomes_selected2gene_expression_list[viral_gn][1])
    for gene_expression in Gene_expression_ME_2020_08_25_10M_B___3300042355_viral_genomes_selected2gene_expression_list[viral_gn][0]:
        if float(gene_expression) > 0:
            total_gene_w_positive_expression_no = total_gene_w_positive_expression_no + 1
    for gene_expression in Gene_expression_ME_2020_08_25_10M_B___3300042355_viral_genomes_selected2gene_expression_list[viral_gn][1]:
        if float(gene_expression) > 0:
            total_gene_w_positive_expression_no = total_gene_w_positive_expression_no + 1
    if float(total_gene_w_positive_expression_no / total_gene_no) < 0.5:
        viral_genomes_not_selected.append(viral_gn)
        
for viral_gn in viral_genomes_not_selected:
    if viral_gn in Gene_expression_ME_2020_08_25_10M_B___3300042355_viral_genomes_selected2gene_expression_list:
        del Gene_expression_ME_2020_08_25_10M_B___3300042355_viral_genomes_selected2gene_expression_list[viral_gn]
        
## Step 4.3 Get the number of viral genomes (with positive expression) and the number of viral genomes that have positive AMG gene expression
no_of_viral_gns_w_positive_expression = 0
no_of_viral_gns_having_positive_AMG_gene_expression = 0
no_of_viral_gns_w_positive_expression = len(Gene_expression_ME_2020_08_25_10M_B___3300042355_viral_genomes_selected2gene_expression_list)
for viral_gn in Gene_expression_ME_2020_08_25_10M_B___3300042355_viral_genomes_selected2gene_expression_list:
    AMG_gene_expression_list = Gene_expression_ME_2020_08_25_10M_B___3300042355_viral_genomes_selected2gene_expression_list[viral_gn][0]
    if any(AMG_gene_expression_list): # The any function returns True if any element of an iterable is True, and False otherwise.
        no_of_viral_gns_having_positive_AMG_gene_expression = no_of_viral_gns_having_positive_AMG_gene_expression + 1

## Step 4.4 Write the result of Step 4.3        
f = open('/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_08_25_10M_B___3300042355_mapping_result_dir/no_of_viral_gn_w_positive_expression_and_AMG_gene_expression.txt','w')
header6 = 'No. of viral genomes with positive expression\tNo. of viral genomes having positive AMG gene expression\n'
f.write(header6)
f.write(f'{no_of_viral_gns_w_positive_expression}\t{no_of_viral_gns_having_positive_AMG_gene_expression}\n')
f.close()

            
# Step 5 Parse the coverage result in 'ME_2020_10_19_5M_B___3300034116_mapping_result_dir' folder   
Coverage_result_ME_2020_10_19_5M_B___3300034116 = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_10_19_5M_B___3300034116_mapping_result_dir/all_coverm_raw_result.txt'
Gene_expression_ME_2020_10_19_5M_B___3300034116_dict = {} # gene => [gene, viral_gene, annotation, gene_expression]
with open(Coverage_result_ME_2020_10_19_5M_B___3300034116, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Ga'):
            tmp = line.split('\t')
            gene, gene_expression = tmp[0], tmp[1]
            if gene in Gene_fun_map: # If the gene is from any viral genomes
                viral_gene, annotation = Gene_fun_map[gene][0], Gene_fun_map[gene][1]
                Gene_expression_ME_2020_10_19_5M_B___3300034116_dict[gene] = [gene, viral_gene, annotation, gene_expression]
lines.close()

## Step 5.1 Write down the gene expression dict
header7 = 'gene\tviral_gene\tannotation\tgene_expression(rpkm)\n'
f = open('/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_10_19_5M_B___3300034116_mapping_result_dir/viral_gene_expression_result.txt','w')
f.write(header7)
for gene in Gene_expression_ME_2020_10_19_5M_B___3300034116_dict:
    line = '\t'.join(Gene_expression_ME_2020_10_19_5M_B___3300034116_dict[gene]) + '\n'
    f.write(line)
f.close()   

## Step 5.2 Get the gene expression result for viral genomes that have at least 50% genes with positive gene expression values
Gene_expression_ME_2020_10_19_5M_B___3300034116_viral_genomes_selected2gene_expression_list = defaultdict(lambda: [[], []]) # viral_gn => [[AMGs],[non-AMGs]]
for gene in Gene_expression_ME_2020_10_19_5M_B___3300034116_dict:   
    viral_gene, annotation, gene_expression = Gene_expression_ME_2020_10_19_5M_B___3300034116_dict[gene][1], Gene_expression_ME_2020_10_19_5M_B___3300034116_dict[gene][2], Gene_expression_ME_2020_10_19_5M_B___3300034116_dict[gene][3]
    viral_gn = viral_gene.split('__')[0] + '__' + viral_gene.split('__')[1]
    if annotation == 'NA':
        Gene_expression_ME_2020_10_19_5M_B___3300034116_viral_genomes_selected2gene_expression_list[viral_gn][1].append(gene_expression)
    elif annotation != 'NA':
        Gene_expression_ME_2020_10_19_5M_B___3300034116_viral_genomes_selected2gene_expression_list[viral_gn][0].append(gene_expression)
        
viral_genomes_not_selected = [] # The viral genomes not pass with positive expression
for viral_gn in Gene_expression_ME_2020_10_19_5M_B___3300034116_viral_genomes_selected2gene_expression_list:
    total_gene_no = 0
    total_gene_w_positive_expression_no = 0
    total_gene_no = len(Gene_expression_ME_2020_10_19_5M_B___3300034116_viral_genomes_selected2gene_expression_list[viral_gn][0]) + len(Gene_expression_ME_2020_10_19_5M_B___3300034116_viral_genomes_selected2gene_expression_list[viral_gn][1])
    for gene_expression in Gene_expression_ME_2020_10_19_5M_B___3300034116_viral_genomes_selected2gene_expression_list[viral_gn][0]:
        if float(gene_expression) > 0:
            total_gene_w_positive_expression_no = total_gene_w_positive_expression_no + 1
    for gene_expression in Gene_expression_ME_2020_10_19_5M_B___3300034116_viral_genomes_selected2gene_expression_list[viral_gn][1]:
        if float(gene_expression) > 0:
            total_gene_w_positive_expression_no = total_gene_w_positive_expression_no + 1
    if float(total_gene_w_positive_expression_no / total_gene_no) < 0.5:
        viral_genomes_not_selected.append(viral_gn)
        
for viral_gn in viral_genomes_not_selected:
    if viral_gn in Gene_expression_ME_2020_10_19_5M_B___3300034116_viral_genomes_selected2gene_expression_list:
        del Gene_expression_ME_2020_10_19_5M_B___3300034116_viral_genomes_selected2gene_expression_list[viral_gn]
        
## Step 5.3 Get the number of viral genomes (with positive expression) and the number of viral genomes that have positive AMG gene expression
no_of_viral_gns_w_positive_expression = 0
no_of_viral_gns_having_positive_AMG_gene_expression = 0
no_of_viral_gns_w_positive_expression = len(Gene_expression_ME_2020_10_19_5M_B___3300034116_viral_genomes_selected2gene_expression_list)
for viral_gn in Gene_expression_ME_2020_10_19_5M_B___3300034116_viral_genomes_selected2gene_expression_list:
    AMG_gene_expression_list = Gene_expression_ME_2020_10_19_5M_B___3300034116_viral_genomes_selected2gene_expression_list[viral_gn][0]
    if any(AMG_gene_expression_list): # The any function returns True if any element of an iterable is True, and False otherwise.
        no_of_viral_gns_having_positive_AMG_gene_expression = no_of_viral_gns_having_positive_AMG_gene_expression + 1

## Step 5.4 Write the result of Step 5.3        
f = open('/storage1/data11/TYMEFLIES_phage/Metatranscriptome/ME_2020_10_19_5M_B___3300034116_mapping_result_dir/no_of_viral_gn_w_positive_expression_and_AMG_gene_expression.txt','w')
header8 = 'No. of viral genomes with positive expression\tNo. of viral genomes having positive AMG gene expression\n'
f.write(header8)
f.write(f'{no_of_viral_gns_w_positive_expression}\t{no_of_viral_gns_having_positive_AMG_gene_expression}\n')
f.close()            