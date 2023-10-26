#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    from collections import defaultdict
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)
    
# Aim: Parse to get the Fst result for summer and winter


def store_seq(file):
    Seq = {}
    head = ""
    with open(file, "r") as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if ' ' in line:
                    head = line.split(' ', 1)[0]
                else:
                    head = line
                Seq[head] = ""
            else:
                Seq[head] += line
    return Seq


# Step 1 Store the Fst result for each viral scf (the cutoff for storing a Fst is 0.15)
viral_scf2fst = {} # viral_scf 2 fst
with open('MetaPop.for_summer_vs_winter/MetaPop/10.Microdiversity/fixation_index.tsv', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('row_samp'):
            tmp = line.split('\t')
            viral_scf, fst = tmp[2], tmp[3]
            if fst != 'NA' and float(fst) > 0.15:
                viral_scf2fst[viral_scf] = fst


# Step 2 Store gene2pi_n_pNpS_info dict
gene2pi_n_pNpS_info = {} # gene => [pi_1, pi_2, pNpS_1, pNpS_2] 
# pi_1 and pi_2 are nucleotide diversity of winter and summer
# pNpS_1 and pNpS_2 are pN/pS ratio of winter and summer
with open('MetaPop.for_summer_vs_winter/MetaPop/10.Microdiversity/global_gene_microdiversity.tsv', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('contig_gene'):
            tmp = line.split('\t')
            gene = tmp[0]  
            gene2pi_n_pNpS_info[gene] = ['NA', 'NA', 'NA', 'NA']
            
with open('MetaPop.for_summer_vs_winter/MetaPop/10.Microdiversity/global_gene_microdiversity.tsv', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('contig_gene'):
            tmp = line.split('\t')
            gene, source, pi, pNpS = tmp[0], tmp[1], tmp[3], tmp[12]            
            if 'winter' in source:
                gene2pi_n_pNpS_info[gene][0] = pi
                if len(pNpS) > 0:
                    gene2pi_n_pNpS_info[gene][2] = pNpS
            elif 'summer' in source:        
                gene2pi_n_pNpS_info[gene][1] = pi
                if len(pNpS) > 0:
                    gene2pi_n_pNpS_info[gene][3] = pNpS                   
                    
                    
# Step 3 Get AMG-containing viral gn list (species representatives)
## Step 3.1 Store species info
Species = {}  # gn_rep => gns
with open("/storage1/data11/TYMEFLIES_phage/Cluster_phage_genomes/Species_level_vOTUs_cluster.txt", "r") as infile:
    for line in infile:
        gn_rep, gns, _ = line.strip().split('\t')
        Species[gn_rep] = gns

## Step 3.2 Store AMG KO information
AMG_summary = {}  # pro => ko
KOs = {}  # ko => ko_detail
IMG2date = {}  # img_id => date_n_season
AMG_containing_viral_gn = {}  # gn => 1
with open("AMG_analysis/AMG_summary.txt", "r") as infile:
    for line in infile:
        if not line.startswith('Pro'):
            tmp = line.strip().split('\t')
            pro, date_n_season, ko, ko_detail = tmp[0], tmp[1], tmp[2], tmp[3]
            AMG_summary[pro] = ko
            img_id = pro.split('__')[0]
            IMG2date[img_id] = date_n_season
            KOs[ko] = ko_detail
            gn = pro.rsplit('__', 1)[0]
            AMG_containing_viral_gn[gn] = 1

### Change the old gene to new gene
Old_gene2new_gene_map = {}  # gene_old => gene_new
New_gene2old_gene_map = {}  # gene_new => gene_old
with open("New_gene2old_gene_map.txt", "r") as infile:
    for line in infile:
        gene_new, gene_old = line.strip().split('\t')
        Old_gene2new_gene_map[gene_old] = gene_new
        New_gene2old_gene_map[gene_new] = gene_old

### Create a list of keys to delete
keys_to_delete = []
for pro in AMG_summary:
    if pro in Old_gene2new_gene_map:
        keys_to_delete.append(pro)

### Iterate over the list of keys to delete and make the changes
for pro in keys_to_delete:
    ko = AMG_summary[pro]
    gene_new = Old_gene2new_gene_map[pro]
    del AMG_summary[pro]  # Delete the old gene and its value
    AMG_summary[gene_new] = ko  # Add the new gene and its value


# Step 4 Get all viral species rep gene annotation
## Step 4.1 Store all the gene ID
All_gene_seq = store_seq("All_phage_species_rep_gn_containing_AMG.mdfed.genes")
All_gene_seq_ID = {}  # Store only the gene ID
for key in All_gene_seq.keys():
    gene = key.replace('>', '', 1)
    All_gene_seq_ID[gene] = 1

All_gene_seq_ID_map = {}  # gene_short_wo_fragment_x => gene_new
for gene in All_gene_seq_ID.keys():
    gene_new = gene

    gene_old = gene
    if gene in New_gene2old_gene_map:
        gene_old = New_gene2old_gene_map[gene]

    gene_short = gene_old.rsplit('__', 1)[1]
    gene_short_wo_fragment_x = gene_short
    if '_fragment_' in gene_short_wo_fragment_x:
        gene_short_wo_fragment_x = gene_short.split('_fragment_', 1)[0] + '_' + gene_short.rsplit('_', 1)[1] # "Ga0334977_0000583_fragment_1_1" will be converted to "Ga0334977_0000583_1"

    All_gene_seq_ID_map[gene_short_wo_fragment_x] = gene_new

## Step 4.2 Store all gene annotation results
annotation_header = ""
All_gene_annotation = {}  # gene => annotation

### Find files matching the pattern
files = glob('/storage1/data11/TYMEFLIES_phage/*/VIBRANT_*.a/VIBRANT_results_*.a/VIBRANT_annotations_*.a.tsv')
for file in files:
    with open(file, 'r') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('protein'):
                annotation_header = '\t'.join(line.split('\t')[1:]) # Delete the first element
            else:
                tmp = line.split('\t')
                gene_short = tmp[0]

                # Get the gene_short_wo_fragment_x
                gene_short_wo_fragment_x = gene_short
                if '_fragment_' in gene_short_wo_fragment_x:
                    gene_short_wo_fragment_x = gene_short.split('_fragment_', 1)[0] + '_' + gene_short.rsplit('_', 1)[1]

                annotation = '\t'.join(line.split('\t')[1:]) # Delete the first element
                if gene_short_wo_fragment_x in All_gene_seq_ID_map:
                    gene = All_gene_seq_ID_map[gene_short_wo_fragment_x]
                    All_gene_annotation[gene] = annotation
                    
                    
# Step 5 Get the final dict of gene2final_fst_results and write it down
gene2final_fst_results = {} # gene => [fst, pi_1, pi_2, pi_compare, pNpS_1, pNpS_2, pNpS_compare] + [annotations] (final fst results; fst indicates the fst of the viral scaffold)
viral_scf2gene_num = defaultdict(int) # viral_scf => gene_num (total gene num of that viral scaffold)
viral_scf2positively_selected_gene_num = defaultdict(int) # viral_scf => positively_selected_gene_num (the number of positively selected genes)
final_fst_result_header_list = ['Gene', 'Fst', 'pi (winter)', 'pi (summer)', 'pi (winter) > pi (summer)', 'pN/pS (winter)', 'pN/pS (summer)', 'pN/pS (winter) < pN/pS (summer)'] + annotation_header.split('\t')
for gene in gene2pi_n_pNpS_info:
    viral_scf = gene.rsplit('_', 1)[0]
    if viral_scf in viral_scf2fst:
        fst = viral_scf2fst[viral_scf]
        pi_1, pi_2, pNpS_1, pNpS_2 = gene2pi_n_pNpS_info[gene]
        pi_compare = 0
        if pi_1 != 'NA' and pi_2 != 'NA' and float(pi_1) > float(pi_2):
            pi_compare = 1
        pNpS_compare = 0
        if pNpS_1 != 'NA' and pNpS_2 != 'NA' and pNpS_1 != 'Inf' and pNpS_2 != 'Inf' and float(pNpS_1) < float(pNpS_2):
            pNpS_compare = 1
        if gene in All_gene_annotation:
            annotations = All_gene_annotation[gene].split('\t')
            gene2final_fst_results[gene] = [fst, pi_1, pi_2, str(pi_compare), pNpS_1, pNpS_2, str(pNpS_compare)] + annotations
        else:
            print(f"{gene} is not in dict of All_gene_annotation")
            
        viral_scf2gene_num[viral_scf] += 1
        if (pi_compare + pNpS_compare) == 2:
            viral_scf2positively_selected_gene_num[viral_scf] += 1            
            
viral_scf2positively_selected_gene_ratio = {} # viral_scf => the percentage of positively selected genes over all genes  
for viral_scf in viral_scf2gene_num:
    gene_num = viral_scf2gene_num[viral_scf] 
    positively_selected_gene_num = viral_scf2positively_selected_gene_num[viral_scf]
    viral_scf2positively_selected_gene_ratio[viral_scf] = positively_selected_gene_num / gene_num
    
## Get the dict for reordering
gene_zfilled2gene = {} # gene_zfilled => gene
for gene in gene2final_fst_results:
    gene_zfilled = gene.rsplit('_', 1)[0] + '_' + gene.rsplit('_', 1)[1].zfill(4)
    gene_zfilled2gene[gene_zfilled] = gene
    
## Write down the gene2final_fst_results 
f = open('MetaPop.for_summer_vs_winter/gene2final_fst_results.txt', 'w')
f.write('\t'.join(final_fst_result_header_list) + '\n')
for gene_zfilled in sorted(gene_zfilled2gene.keys()):
    gene = gene_zfilled2gene[gene_zfilled]
    final_fst_results = gene2final_fst_results[gene]
    f.write(gene + '\t' + '\t'.join(final_fst_results) + '\n')
f.close()    

## Write down the viral_scf2positively_selected_gene_ratio
f = open('MetaPop.for_summer_vs_winter/viral_scf2positively_selected_gene_ratio.txt', 'w')
f.write('viral scaffold\tpositively selected gene ratio\n')
for viral_scf in sorted(viral_scf2positively_selected_gene_ratio.keys()):
    line = viral_scf + '\t' + str(viral_scf2positively_selected_gene_ratio[viral_scf])
    f.write(line + '\n')
f.close()    
