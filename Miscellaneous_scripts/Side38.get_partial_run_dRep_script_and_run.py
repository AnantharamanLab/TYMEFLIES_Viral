#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    from collections import defaultdict
    from pathlib import Path
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: Get partial (left) run_dRep script from the list of dRep output folders


# Step 1 Store the list of dRep output folders
dRep_output_folder_list = []
with open('tmp.dRep_output_folder_list.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        dRep_output_folder_list.append(line)

# Step 2 Store the partial run_dRep script
run_dRep_script_list = []
with open('tmp.run_dRep.sh', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.split(' ')[2].split('/')[1] not in dRep_output_folder_list:        
            run_dRep_script_list.append(line)


# Step 3 Write down the partial run_dRep script and run it 
f = open('tmp.run_dRep.partial.sh', 'w')
for line in run_dRep_script_list:
    f.write(line + '\n')
f.close()

os.system("cat tmp.run_dRep.partial.sh | parallel -j 20")
os.system("rm tmp.run_dRep.partial.sh")  

       