#!/usr/bin/env python3

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
    
# Aim: Filter AMGs in an integrative manner and add metagenome and KEGG info to AMG summary table
# There are following steps to filter AMGs
# (1) Filter tail AMGs
# (2) Filter AMGs that have any v-scores (KEGG and Pfam v-scores) >= 1
# (3) Filter AMGs with flanking genes of KEGG v-scores < 0.25
# (4) Filter AMGs that have COG category as T or B


def store_seq(input_seq_file): # The input sequence file should be a file with full path
    head = "" # Store the header line
    seq_dict = {} # Store the sequence dict
    
    with open(input_seq_file, "r") as seq_lines:
        for line in seq_lines:
            line = line.rstrip("\n") # Remove "\n" in the end
            if ">" in line:
                head = line.split(None, 1)[0] # Cut at the first " " or "\t", use the first part
                seq_dict[head] = ""                 
            else:
                seq_dict[head] += line
            
    seq_lines.close()
    
    return seq_dict
    
def check_left(AMG, pros_ordered, AMG_list):
    AMG_index = pros_ordered.index(AMG) # The index of the given AMG
    logic = True
    for i in range(0, AMG_index):
        if pros_ordered[i] not in AMG_list:
            logic = False
    return logic   

def check_right(AMG, pros_ordered, AMG_list):
    AMG_index = pros_ordered.index(AMG) # The index of the given AMG
    logic = True
    for i in range((AMG_index + 1), len(pros_ordered)):
        if pros_ordered[i] not in AMG_list:
            logic = False
    return logic  
    
    
# Step 1 Filter tail AMGs
## Step 1.1 Open the file containing the AMG list and read the AMGs into a list
All_AMG_list = [] # Store all the AMG list
with open('AMG_analysis/AMG_summary.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('Protein'):
            AMG = line.split('\t')[0]
            All_AMG_list.append(AMG)

### Store the scf2AMG_list dict
scf2AMG_list = defaultdict(list)
for AMG in All_AMG_list:
    scf = AMG.rsplit('_', 1)[0]
    scf2AMG_list[scf].append(AMG)

## Step 1.2 Store scaffold to protein list dict
scf2pros_ordered = {} # scf => [pros]; The proteins are re-ordered
scf2pros = defaultdict(list) # scf => [pros]

all_virus_protein_seq = {} # The seq dict for all virus proteins
all_viral_faa_addrs = glob('/storage1/data11/TYMEFLIES_phage/*/vRhyme_best_bins_fasta_parsed/*.faa')
for each_viral_faa_addr in all_viral_faa_addrs:
    each_viral_faa_seq = store_seq(each_viral_faa_addr)
    all_virus_protein_seq.update(each_viral_faa_seq)

headers = [header.replace('>', '', 1) for header in all_virus_protein_seq]
for header in headers:
    scf = header.rsplit('_', 1)[0]
    pro = header
    scf2pros[scf].append(pro)
    
for scf in scf2pros:
    pros = scf2pros[scf] # pros is now a list
    # Make dict from [pros] as: pro => order_num
    pro2order_num = {pro:int(pro.rsplit('_', 1)[1]) for pro in pros}
    
    pros_ordered = list(sorted(pro2order_num, key = pro2order_num.get)) # pros_ordered is now a re-ordered list
    scf2pros_ordered[scf] = pros_ordered
    
## Step 1.3 Get the tail AMG label results
AMG2label = {} # AMG => [position, tail_AMG (or not_tail_AMG)]; position can be '1 in 11'
for scf in scf2AMG_list:
    pros_ordered = scf2pros_ordered[scf]
    AMG_list = scf2AMG_list[scf]
    
    # First, start from left to right
    for AMG in AMG_list:
        AMG_index = pros_ordered.index(AMG) # The index of the given AMG
        if AMG_index == 0:
            AMG2label[AMG] = [f'{(AMG_index + 1)} in {len(pros_ordered)}', 'tail_AMG']
        if AMG_index > 0:
            if check_left(AMG, pros_ordered, AMG_list): # If the left genes are all AMGs 
                AMG2label[AMG] = [f'{(AMG_index + 1)} in {len(pros_ordered)}', 'tail_AMG']
     
    # Second, start from right to left
    for AMG in AMG_list:
        AMG_index = pros_ordered.index(AMG) # The index of the given AMG
        if AMG_index == (len(pros_ordered) - 1) :
            AMG2label[AMG] = [f'{(AMG_index + 1)} in {len(pros_ordered)}', 'tail_AMG']
        if AMG_index < (len(pros_ordered) - 1):
            if check_right(AMG, pros_ordered, AMG_list): # If the left genes are all AMGs 
                AMG2label[AMG] = [f'{(AMG_index + 1)} in {len(pros_ordered)}', 'tail_AMG']
    
    # Label the rest as 'not_tail_AMG'
    for AMG in AMG_list:
        AMG_index = pros_ordered.index(AMG) # The index of the given AMG
        if AMG not in AMG2label:
            AMG2label[AMG] = [f'{(AMG_index + 1)} in {len(pros_ordered)}', 'not_tail_AMG']
                
### Filter tail_AMG
tail_AMG = set()
for AMG in AMG2label:
    if AMG2label[AMG][1] == 'tail_AMG':
        tail_AMG.add(AMG)
        
       
# Step 2 Filter AMGs that have any v-scores (KEGG and Pfam v-scores) >= 1   
## Step 2.1 Store all viral protein annotation result
VIBRANT_annotation_result_addrs = glob('/storage1/data11/TYMEFLIES_phage/*/VIBRANT_*.a/VIBRANT_results_*.a/VIBRANT_annotations_*.a.tsv')
pro_short2VIBRANT_annotation = {} # pro_short => [VIBRANT_annotation]
for VIBRANT_annotation_result_addr in VIBRANT_annotation_result_addrs:
    with open(VIBRANT_annotation_result_addr, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if line.startswith('Protein'):
                header = line
            else:
                tmp = line.split('\t')
                pro_short = tmp[0]
                pro_short2VIBRANT_annotation[pro_short] = tmp
    lines.close()
    
### Store pro_short_w_fragment2VIBRANT_annotation dict
pro_short_w_fragment2VIBRANT_annotation = {} # pro_short_w_fragment => [VIBRANT_annotation]
for pro_short in pro_short2VIBRANT_annotation:
    if '_fragment_' in pro_short:   
        pro_short_w_fragment2VIBRANT_annotation[pro_short] = pro_short2VIBRANT_annotation[pro_short]
       
## Step 2.2 Replace short protein to long protein for all viral proteins    
Pro2VIBRANT_annotation = {} # pro => [VIBRANT_annotation]    
for pro in headers:   
    pro_short_from_pro = pro.rsplit('__', 1)[1]
    if pro_short_from_pro in pro_short2VIBRANT_annotation:
        Pro2VIBRANT_annotation[pro] = pro_short2VIBRANT_annotation[pro_short_from_pro]
    else: # Need to delete fragment within the pro name
        for pro_short in pro_short_w_fragment2VIBRANT_annotation: # All the pro_short in this dict contain '_fragment_' inside 
                pro_short_mdfed = pro_short.split('_fragment_', 1)[0] + '_' + pro_short.rsplit('_', 1)[1] # Delete '_fragment_{i}_' from the pro_short
                if pro_short_mdfed == pro_short_from_pro:
                    Pro2VIBRANT_annotation[pro] = pro_short_w_fragment2VIBRANT_annotation[pro_short]
                   
## Step 2.3 Filter AMGs based on v-scores
AMGs_w_high_v_scores = [] # Store the list of AMGs that have any v-scores (KEGG and Pfam v-scores) >= 1
for AMG in All_AMG_list:
    KO_v_score, Pfam_v_score = Pro2VIBRANT_annotation[AMG][7], Pro2VIBRANT_annotation[AMG][12]
    logic = True
    if KO_v_score and float(KO_v_score) >= 1:
        logic = False
    if Pfam_v_score and float(Pfam_v_score) >= 1:
        logic = False        
    if not logic:
        AMGs_w_high_v_scores.append(AMG)


# Step 3 Filter AMGs with flanking genes of v-scores < 0.25
## Step 3.1 Get AMG to flanking genes dict
flanking_num = 2 # The number of flanking gene on both sides to be taken into consideration
AMG2flanking_genes = {} # AMG => [flanking genes]
        
for scf in scf2AMG_list:
    pros_ordered = scf2pros_ordered[scf]
    AMG_list = scf2AMG_list[scf]
    
    for AMG in AMG_list:
        if AMG2label[AMG][1] == 'not_tail_AMG': # Only focusing not_tail_AMG
            flanking_genes_left = []
            flanking_genes_right = []
            AMG_index = pros_ordered.index(AMG) # The index of the given AMG
            # Find left side flanking genes
            for i in range((AMG_index - 1), -1, -1):
                if pros_ordered[i] not in AMG_list:
                    flanking_genes_left.append(pros_ordered[i])
                    if len(flanking_genes_left) == flanking_num:
                        break
            # Find right side flanking genes:
            for i in range((AMG_index + 1), len(pros_ordered), 1):
                if pros_ordered[i] not in AMG_list:
                    flanking_genes_right.append(pros_ordered[i])
                    if len(flanking_genes_right) == flanking_num:
                        break  
            flanking_genes = flanking_genes_left + flanking_genes_right
            AMG2flanking_genes[AMG] = flanking_genes
            
## Step 3.2 Filter AMGs with flanking genes of v-scores < 0.25
AMGs_w_flanking_genes_of_low_v_score = []
for AMG in AMG2flanking_genes:
    flanking_genes = AMG2flanking_genes[AMG]
    low_v_score_gene_count = 0
    for gene in flanking_genes:
        if gene in Pro2VIBRANT_annotation:
            KO_v_score = Pro2VIBRANT_annotation[gene][7]
            if (KO_v_score and float(KO_v_score) < 0.25):
                low_v_score_gene_count = low_v_score_gene_count + 1 

    if low_v_score_gene_count == len(flanking_genes): # If all the flanking genes are low v-score genes
        AMGs_w_flanking_genes_of_low_v_score.append(AMG) 
        
   
# Step 4 Filter AMGs that have COG category as T or B
## Step 4.1 Store the eggNOG result
AMG2eggNOG_annotation = {} # AMG => [eggNOG_annotation] 
with open('AMG.emapper.annotations.tsv', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('#'):
            tmp = line.split('\t')
            AMG2eggNOG_annotation[tmp[0]] = tmp
lines.close()

## Step 4.2 Filter AMGs based on COG category
AMGs_belong_to_not_correct_COG = [] # Store the list of AMGs that have COG category as T or B
for AMG in AMG2eggNOG_annotation:
    eggNOG_annotation = AMG2eggNOG_annotation[AMG]
    COG_category = eggNOG_annotation[6]
    if 'T' in COG_category or 'B' in COG_category:
        AMGs_belong_to_not_correct_COG.append(AMG)


# Step 5 Write down the final filter AMGs
AMG_filtered = set()
AMG_filtered = set(All_AMG_list) - tail_AMG - set(AMGs_w_high_v_scores) - set(AMGs_w_flanking_genes_of_low_v_score) - set(AMGs_belong_to_not_correct_COG)

print (f"total AMG number before filtering is {len(All_AMG_list)}")
print (f"total AMG number after filtering is {len(AMG_filtered)}")
print (f"tail_AMG number is {len(tail_AMG)}")
print (f"AMGs_w_high_v_scores number is {len(AMGs_w_high_v_scores)}")
print (f"AMGs_w_flanking_genes_of_low_v_score number is {len(AMGs_w_flanking_genes_of_low_v_score)}")
print (f"AMGs_belong_to_not_correct_COG number is {len(AMGs_belong_to_not_correct_COG)}")
                                  
f = open('AMG_analysis/AMG_filtered_list.txt', 'w')
for AMG in AMG_filtered:
    f.write(AMG + '\n')
f.close()


# Step 6 Store all metagenomes info
Meta_info = {}  # img_id => date_and_season (for example, "2000-03-15	| 75 | Spring")
with open("TYMEFLIES_metagenome_info.txt") as f:
    for line in f:
        line = line.rstrip('\n')
        if not line.startswith("IMG"):
            tmp = line.split("\t")
            img_id = tmp[0]
            date_and_season = f"{tmp[8]} | {tmp[9]} | {tmp[10]}"
            Meta_info[img_id] = date_and_season

KEGG_pathway_info = {}  # map_ => [0] metabolism [1] pathway [2] kos
with open("VIBRANT_KEGG_pathways_summary.tsv") as f:
    for line in f:
        line = line.rstrip('\n')
        if not line.startswith("Entry"):
            tmp = line.split("\t")
            map_ = tmp[0]
            metabolism = tmp[1]
            pathway = tmp[2]
            kos = tmp[3]
            KEGG_pathway_info[map_] = [metabolism, pathway, kos]

KEGG_module_info = {}  # module_step => [0] module [1] kos
with open("/slowdata/data1/Genome_profile_software/METABOLIC_template_and_database/kegg_module_step_db.txt") as f:
    for line in f:
        line = line.rstrip('\n')
        if not line.startswith("name"):
            tmp = line.split("\t")
            module = tmp[0]
            kos = tmp[1]
            module_step = tmp[2]
            KEGG_module_info[module_step] = [module, kos]

AMG_summary = {}  # pro => [0] amg_ko [1] amg_ko_name [2] pfam [3] pfam_name
with open("AMG_analysis/AMG_summary.txt") as f:
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('Protein'):
            # skip header line
            continue
        else:
            tmp = line.split("\t")
            pro = tmp[0]
            amg_ko = tmp[1]
            amg_ko_name = tmp[2]
            if amg_ko_name.startswith('"') and amg_ko_name.endswith('"'):
                # Remove the quote marks
                amg_ko_name = amg_ko_name[1:-1]
            pfam = tmp[3]
            pfam_name = tmp[4]
            if pfam_name.startswith('"') and pfam_name.endswith('"'):
                # Remove the quote marks
                pfam_name = pfam_name[1:-1]            
            AMG_summary[pro] = [amg_ko, amg_ko_name, pfam, pfam_name]
            
## Add metagenome info, KEGG pathway info, KEGG module info to the AMG summary dict
for pro in sorted(AMG_summary.keys()):
    amg_ko = AMG_summary[pro][0]
    amg_ko_name = AMG_summary[pro][1]
    pfam = AMG_summary[pro][2]
    pfam_name = AMG_summary[pro][3]

    img_id = pro.split("__")[0]
    date_and_season = Meta_info[img_id]

    metabolisms = ""
    pathways = ""
    modules = ""

    Metabolisms = {}
    Pathways = {}
    Modules = {}

    for map_ in sorted(KEGG_pathway_info.keys()):
        kos = KEGG_pathway_info[map_][2]
        if amg_ko in kos:
            metabolism = KEGG_pathway_info[map_][0]
            pathway = KEGG_pathway_info[map_][1]

            Metabolisms[metabolism] = 1
            Pathways[pathway] = 1

    for module_step in sorted(KEGG_module_info.keys()):
        kos = KEGG_module_info[module_step][1]
        if amg_ko in kos:
            module = KEGG_module_info[module_step][0]

            Modules[module] = 1

    if len(Metabolisms) > 0:
        metabolisms = "|".join(sorted(Metabolisms.keys()))
    if len(Pathways) > 0:
        pathways = "|".join(sorted(Pathways.keys()))
    if len(Modules) > 0:
        modules = "|".join(sorted(Modules.keys()))

    AMG_summary[pro] = AMG_summary[pro] + [date_and_season, metabolisms, pathways, modules]


## Write down new AMG_summary.txt
with open("AMG_analysis/AMG_summary_new.txt", "w") as outfile:
    outfile.write("Protein\tMetagenome date and season ('date | date in the year | season')\tAMG KO\tAMG KO name\tPfam\tPfam name\tMetabolisms (mutiple metabolisms separated by '|')\tPathways (mutiple pathways separated by '|')\tModules (mutiple modules separated by '|')\n")
    for pro in AMG_filtered:
        date_and_season = AMG_summary[pro][4]
        amg_ko = AMG_summary[pro][0]
        amg_ko_name = AMG_summary[pro][1]
        pfam = AMG_summary[pro][2]
        pfam_name = AMG_summary[pro][3]
        metabolisms = AMG_summary[pro][5]
        pathways = AMG_summary[pro][6]
        modules = AMG_summary[pro][7]
        outfile.write(f"{pro}\t{date_and_season}\t{amg_ko}\t{amg_ko_name}\t{pfam}\t{pfam_name}\t{metabolisms}\t{pathways}\t{modules}\n")
            
## Sort the result by date
os.system('cat AMG_analysis/AMG_summary_new.txt | sort -k 2 -n > tmp')
os.system('mv tmp AMG_analysis/AMG_summary_new.txt')

# Replace the old summary file
os.system('rm AMG_analysis/AMG_summary.txt')
os.system('mv AMG_analysis/AMG_summary_new.txt AMG_analysis/AMG_summary.txt')           
