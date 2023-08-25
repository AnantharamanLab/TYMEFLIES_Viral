#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
# Aim: Make the mapping reference for 2020 metatranscriptomic datasets 
# The metaT ID to corresponding metaG ID:
# ME_2020_07_24_5m_C => 3300044631 (2018/7/24)
# ME_2020_08_05_10M_B => 3300044608 (2018/8/6)
# ME_2020_08_25_10M_B => 3300042355 (2018/8/25)
# ME_2020_10_19_5M_B => 3300034116 (2018/10/24)

# Use 3300044631 (2018/7/24) as the mapping reference for ME_2020_07_24_5m_C 
# Use 3300044608 (2018/8/6) as the mapping reference for ME_2020_08_05_10M_B 
# Use 3300042355 (2018/8/25) as the mapping reference for ME_2020_08_25_10M_B 
# Use 3300034116 (2018/10/24) as the mapping reference for ME_2020_10_19_5M_B

def store_seq(input_seq_file): # The input sequence file should be a file with full path
    head = '' # Store the header line
    seq_dict = {} # Store the sequence dict
    
    with open(input_seq_file, 'r') as seq_lines:
        for line in seq_lines:
            line = line.rstrip('\n') # Remove '\n' in the end
            if '>' in line:
                if (' ' or '\t') in line: # Break at the first ' ' or '\t'
                    spliter = ''
                    for i in range(len(line)):
                        if line[i] == ' ' or line[i] == '\t':
                            spliter = line[i]
                            break 
                           
                    head = line.split(f'{spliter}', 1)[0]
                    seq_dict[head] = ''
                else:
                    head = line
                    seq_dict[head] = ''
            else:
                seq_dict[head] += line
            
    seq_lines.close()
    
    return seq_dict


# Step 1 Store the AMG map
AMG_map = {} # amg_gene => [KO, AMG_KO_name]; Store the annotation of all AMG genes
with open('/storage1/data11/TYMEFLIES_phage/AMG_analysis/AMG_summary.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('Protein'):
            tmp = line.split('\t')
            amg_gene, KO, AMG_KO_name = tmp[0], tmp[2], tmp[3]
            AMG_map[amg_gene] = [KO, AMG_KO_name]
lines.close()            


Gene_fun_map = {} # gene => [(NA or viral_gene_id), annotation]; Store maps for both metagenomes

# Step 2 Get 3300044631 all genes and all gene function map
IMG3300044631_ffn_seq = store_seq('/storage1/data11/TYMEFLIES_phage/3300044631/VIBRANT_3300044631.a/3300044631.a.prodigal.ffn')

IMG3300044631_all_vRhyme_bin_ffn_addrs = glob('/storage1/data11/TYMEFLIES_phage/3300044631/vRhyme_best_bins_fasta_parsed/*.ffn')
IMG3300044631_all_vRhyme_bin_gene2gene = {} # Store viral bin gene ID to original gene ID map
Gene2IMG3300044631_all_vRhyme_bin_gene = {} # Store original gene ID to viral bin gene ID map
for each_vRhyme_bin_ffn_addr in IMG3300044631_all_vRhyme_bin_ffn_addrs:
    each_vRhyme_bin_ffn_seq = store_seq(each_vRhyme_bin_ffn_addr)
    for header in each_vRhyme_bin_ffn_seq:
        viral_bin_gene_ID = header.replace('>', '', 1)
        original_ID = viral_bin_gene_ID.rsplit('__', 1)[1]
        if '_fragment_' in original_ID:
            original_ID = original_ID.split('_fragment', 1)[0] + original_ID.rsplit('_', 1)[1]
        IMG3300044631_all_vRhyme_bin_gene2gene[viral_bin_gene_ID] = original_ID   
        Gene2IMG3300044631_all_vRhyme_bin_gene[original_ID] = viral_bin_gene_ID        
            
for header in IMG3300044631_ffn_seq:
    gene = header.replace('>', '', 1)
    viral_gene_id = 'NA'
    if gene in Gene2IMG3300044631_all_vRhyme_bin_gene:
        viral_gene_id = Gene2IMG3300044631_all_vRhyme_bin_gene[gene]
    annotation = 'NA'
    if viral_gene_id in AMG_map:
        annotation = AMG_map[viral_gene_id][0] + ' | ' + AMG_map[viral_gene_id][1]
    Gene_fun_map[gene] = [viral_gene_id, annotation]


# Step 3 Get 3300044608 all genes and all gene function map
IMG3300044608_ffn_seq = store_seq('/storage1/data11/TYMEFLIES_phage/3300044608/VIBRANT_3300044608.a/3300044608.a.prodigal.ffn')

IMG3300044608_all_vRhyme_bin_ffn_addrs = glob('/storage1/data11/TYMEFLIES_phage/3300044608/vRhyme_best_bins_fasta_parsed/*.ffn')
IMG3300044608_all_vRhyme_bin_gene2gene = {} # Store viral bin gene ID to original gene ID map
Gene2IMG3300044608_all_vRhyme_bin_gene = {} # Store original gene ID to viral bin gene ID map
for each_vRhyme_bin_ffn_addr in IMG3300044608_all_vRhyme_bin_ffn_addrs:
    each_vRhyme_bin_ffn_seq = store_seq(each_vRhyme_bin_ffn_addr)
    for header in each_vRhyme_bin_ffn_seq:
        viral_bin_gene_ID = header.replace('>', '', 1)
        original_ID = viral_bin_gene_ID.rsplit('__', 1)[1]
        if '_fragment_' in original_ID:
            original_ID = original_ID.split('_fragment', 1)[0] + original_ID.rsplit('_', 1)[1]
        IMG3300044608_all_vRhyme_bin_gene2gene[viral_bin_gene_ID] = original_ID   
        Gene2IMG3300044608_all_vRhyme_bin_gene[original_ID] = viral_bin_gene_ID        
            
for header in IMG3300044608_ffn_seq:
    gene = header.replace('>', '', 1)
    viral_gene_id = 'NA'
    if gene in Gene2IMG3300044608_all_vRhyme_bin_gene:
        viral_gene_id = Gene2IMG3300044608_all_vRhyme_bin_gene[gene]
    annotation = 'NA'
    if viral_gene_id in AMG_map:
        annotation = AMG_map[viral_gene_id][0] + ' | ' + AMG_map[viral_gene_id][1]
    Gene_fun_map[gene] = [viral_gene_id, annotation]


# Step 4 Get 3300042355 all genes and all gene function map
IMG3300042355_ffn_seq = store_seq('/storage1/data11/TYMEFLIES_phage/3300042355/VIBRANT_3300042355.a/3300042355.a.prodigal.ffn')

IMG3300042355_all_vRhyme_bin_ffn_addrs = glob('/storage1/data11/TYMEFLIES_phage/3300042355/vRhyme_best_bins_fasta_parsed/*.ffn')
IMG3300042355_all_vRhyme_bin_gene2gene = {} # Store viral bin gene ID to original gene ID map
Gene2IMG3300042355_all_vRhyme_bin_gene = {} # Store original gene ID to viral bin gene ID map
for each_vRhyme_bin_ffn_addr in IMG3300042355_all_vRhyme_bin_ffn_addrs:
    each_vRhyme_bin_ffn_seq = store_seq(each_vRhyme_bin_ffn_addr)
    for header in each_vRhyme_bin_ffn_seq:
        viral_bin_gene_ID = header.replace('>', '', 1)
        original_ID = viral_bin_gene_ID.rsplit('__', 1)[1]
        if '_fragment_' in original_ID:
            original_ID = original_ID.split('_fragment', 1)[0] + original_ID.rsplit('_', 1)[1]
        IMG3300042355_all_vRhyme_bin_gene2gene[viral_bin_gene_ID] = original_ID   
        Gene2IMG3300042355_all_vRhyme_bin_gene[original_ID] = viral_bin_gene_ID        
            
for header in IMG3300042355_ffn_seq:
    gene = header.replace('>', '', 1)
    viral_gene_id = 'NA'
    if gene in Gene2IMG3300042355_all_vRhyme_bin_gene:
        viral_gene_id = Gene2IMG3300042355_all_vRhyme_bin_gene[gene]
    annotation = 'NA'
    if viral_gene_id in AMG_map:
        annotation = AMG_map[viral_gene_id][0] + ' | ' + AMG_map[viral_gene_id][1]
    Gene_fun_map[gene] = [viral_gene_id, annotation]
    
    
# Step 5 Get 3300034116 all genes and all gene function map
IMG3300034116_ffn_seq = store_seq('/storage1/data11/TYMEFLIES_phage/3300034116/VIBRANT_3300034116.a/3300034116.a.prodigal.ffn')

IMG3300034116_all_vRhyme_bin_ffn_addrs = glob('/storage1/data11/TYMEFLIES_phage/3300034116/vRhyme_best_bins_fasta_parsed/*.ffn')
IMG3300034116_all_vRhyme_bin_gene2gene = {} # Store viral bin gene ID to original gene ID map
Gene2IMG3300034116_all_vRhyme_bin_gene = {} # Store original gene ID to viral bin gene ID map
for each_vRhyme_bin_ffn_addr in IMG3300034116_all_vRhyme_bin_ffn_addrs:
    each_vRhyme_bin_ffn_seq = store_seq(each_vRhyme_bin_ffn_addr)
    for header in each_vRhyme_bin_ffn_seq:
        viral_bin_gene_ID = header.replace('>', '', 1)
        original_ID = viral_bin_gene_ID.rsplit('__', 1)[1]
        if '_fragment_' in original_ID:
            original_ID = original_ID.split('_fragment', 1)[0] + '_' + original_ID.rsplit('_', 1)[1]
        IMG3300034116_all_vRhyme_bin_gene2gene[viral_bin_gene_ID] = original_ID   
        Gene2IMG3300034116_all_vRhyme_bin_gene[original_ID] = viral_bin_gene_ID        
            
for header in IMG3300034116_ffn_seq:
    gene = header.replace('>', '', 1)
    viral_gene_id = 'NA'
    if gene in Gene2IMG3300034116_all_vRhyme_bin_gene:
        viral_gene_id = Gene2IMG3300034116_all_vRhyme_bin_gene[gene]
    annotation = 'NA'
    if viral_gene_id in AMG_map:
        annotation = AMG_map[viral_gene_id][0] + ' | ' + AMG_map[viral_gene_id][1]
    Gene_fun_map[gene] = [viral_gene_id, annotation]

    
# Step 6 Concatenate to make the final mapping reference and make the final gene function map
map_ref_file1 = '/storage1/data11/TYMEFLIES_phage/3300044631/VIBRANT_3300044631.a/3300044631.a.prodigal.ffn'
map_ref_file2 = '/storage1/data11/TYMEFLIES_phage/3300044608/VIBRANT_3300044608.a/3300044608.a.prodigal.ffn'
map_ref_file3 = '/storage1/data11/TYMEFLIES_phage/3300042355/VIBRANT_3300042355.a/3300042355.a.prodigal.ffn'
map_ref_file4 = '/storage1/data11/TYMEFLIES_phage/3300034116/VIBRANT_3300034116.a/3300034116.a.prodigal.ffn'
the_2020_metaT_ref_fa_for_3300044631 = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_2020_metaT_ref_for_3300044631.based_on_MetaG.fa'
the_2020_metaT_ref_fa_for_3300044608 = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_2020_metaT_ref_for_3300044608.based_on_MetaG.fa'
the_2020_metaT_ref_fa_for_3300042355 = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_2020_metaT_ref_for_3300042355.based_on_MetaG.fa'
the_2020_metaT_ref_fa_for_3300034116 = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_2020_metaT_ref_for_3300034116.based_on_MetaG.fa'
os.system(f'cat {map_ref_file1} > {the_2020_metaT_ref_fa_for_3300044631}')
os.system(f'cat {map_ref_file2} > {the_2020_metaT_ref_fa_for_3300044608}')
os.system(f'cat {map_ref_file3} > {the_2020_metaT_ref_fa_for_3300042355}')
os.system(f'cat {map_ref_file4} > {the_2020_metaT_ref_fa_for_3300034116}')

the_2020_metaT_ref_gene_function_map = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_2020_metaT_ref_gene_function_map.based_on_MetaG.txt'
f = open(the_2020_metaT_ref_gene_function_map,'w')
for gene in sorted(Gene_fun_map.keys()):
    line = gene + '\t' + Gene_fun_map[gene][0] + '\t' + Gene_fun_map[gene][1] + '\n'
    f.write(line)
f.close()    


