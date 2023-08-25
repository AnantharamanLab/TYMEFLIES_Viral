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
    
# Aim: Make the mapping reference for 24 hrs metatranscriptomic datasets (2015/8/20-21)
#      Use All_phage_species_rep_gn_containing_AMG.genes and AMG_counterpart_genes_and_flankings.fasta as the mapping reference

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


# Step 2 Concatenate to make the final mapping reference and make the final gene function map
All_phage_species_rep_gn_gene_fun_map = {} # gene => NA or 'KO | AMG_KO_name'
map_ref_file1 = '/storage1/data11/TYMEFLIES_phage/All_phage_species_rep_gn_containing_AMG.genes'
map_ref_file2 = '/storage1/data11/TYMEFLIES_phage/AMG_counterpart_genes_and_flankings.fasta'
the_24hrs_metaT_ref_fa = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_24hrs_metaT_ref.based_on_all_phage_species_rep.fa'
os.system(f'cat {map_ref_file1} {map_ref_file2} > {the_24hrs_metaT_ref_fa}')

the_24hrs_metaT_ref_gene_function_map = '/storage1/data11/TYMEFLIES_phage/Metatranscriptome/the_24hrs_metaT_ref_gene_function_map.based_on_all_phage_species_rep.txt'
All_phage_species_rep_gn_gene_seq = store_seq(map_ref_file1)
for header in All_phage_species_rep_gn_gene_seq:
    gene = header.replace('>', '', 1)
    if gene in AMG_map:
        All_phage_species_rep_gn_gene_fun_map[gene] = AMG_map[gene][0] + ' | ' + AMG_map[gene][1]
    else:
        All_phage_species_rep_gn_gene_fun_map[gene] = 'NA'

f = open(the_24hrs_metaT_ref_gene_function_map,'w')
for gene in sorted(All_phage_species_rep_gn_gene_fun_map.keys()):
    line = gene + '\t' + All_phage_species_rep_gn_gene_fun_map[gene] + '\n'
    f.write(line)
f.close()    