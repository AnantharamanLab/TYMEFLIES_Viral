#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    warnings.filterwarnings("ignore")
    from pathlib import Path
    from collections import defaultdict    
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)

# Aim: make vOTU representative list; pick the high quality and long virus genome

# Step 1 Store all genome property
Gn_property = {} # $gn => [vOTU, length, completeness, genome_quality, tax]
vOTU2gns = defaultdict(list) # $vOTU => [gns]
with open("IMGVR_all_Sequence_information.tsv", 'r') as IN:
    for line in IN:
        if not line.startswith('UVIG'):
            tmp = line.rstrip('\n').split('\t')
            gn = tmp[0]
            vOTU = tmp[5]
            length = tmp[6]
            completeness = tmp[10]
            genome_quality = tmp[12]
            tax = tmp[14]
            Gn_property[gn] = [vOTU, length, completeness, genome_quality, tax]
            vOTU2gns[vOTU].append(gn) # Add gn list in each vOTU
IN.close()            

# Step 2 Pick vOTU representative 
vOTU2rep = {} # vOTU => gn_rep (representative genome)
for vOTU in sorted(vOTU2gns.keys()):
    Gns = vOTU2gns[vOTU] # The list that contains all gns in this vOTU
    
    # Divide genomes into 4 collections by genome quality
    Gns2 = [] # 2nd level genome collection for picking, containing "Reference"
    Gns3 = [] # 3rd level genome collection for picking, containing "High-quality"
    Gns4 = [] # 4th level genome collection for picking, containing "Genome fragment"
    Gns5 = [] # 5th level genome collection for picking, containing "Unsure" (>120% or no completeness estimate) 
    for gn in Gns:
        if Gn_property[gn][3] == "Reference":
            Gns2.append(gn)
        elif Gn_property[gn][3] == "High-quality":
            Gns3.append(gn)
        elif Gn_property[gn][3] == "Genome fragment":
            Gns4.append(gn)
        elif Gn_property[gn][3] == "Unsure":
            Gns5.append(gn)

    if Gns2: # If the 2nd level genome collection is not empty
        gn_rep = Gns2[0]
        for i in range(1, len(Gns2)):
            if int(Gn_property[Gns2[i]][1]) > int(Gn_property[gn_rep][1]): # Pick the longest one from the 2nd level genome collection
                gn_rep = Gns2[i]
        vOTU2rep[vOTU] = gn_rep
    elif not Gns2 and Gns3: # If the 2nd level genome collection is empty and the 3rd level genome collection is not empty
        gn_rep = Gns3[0]
        for i in range(1, len(Gns3)):
            if int(Gn_property[Gns3[i]][1]) > int(Gn_property[gn_rep][1]): # Pick the longest one from the 3rd level genome collection
                gn_rep = Gns3[i]
        vOTU2rep[vOTU] = gn_rep
    elif not Gns2 and not Gns3 and Gns4: # If the 2nd and 3rd level genome collection are empty and the 4th level genome collection is not empty
        gn_rep = Gns4[0]
        for i in range(1, len(Gns4)):
            if int(Gn_property[Gns4[i]][1]) > int(Gn_property[gn_rep][1]): # Pick the longest one from the 4th level genome collection
                gn_rep = Gns4[i]
        vOTU2rep[vOTU] = gn_rep
    elif not Gns2 and not Gns3 and not Gns4 and Gns5: # If all the above genome collections are empty, pick the first one from the 5th level genome collection
        vOTU2rep[vOTU] = Gns5[0]
    else:
        warnings.warn("No genomes found for vOTU: " + vOTU)

# Step 3 Write the vOTU representative list
with open("IMGVR_all_phage_vOTU_representatives.txt", 'w') as OUT:
    OUT.write("vOTU\tvOTU representative\tvOTU representative taxonomy\n")
    for vOTU, gn_rep in vOTU2rep.items():
        OUT.write(vOTU + "\t" + gn_rep + "\t" + Gn_property[gn_rep][4] + "\n")
OUT.close()  

