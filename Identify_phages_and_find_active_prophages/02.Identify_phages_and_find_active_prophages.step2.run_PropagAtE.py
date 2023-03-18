#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    from pathlib import Path 
    from collections import defaultdict
    warnings.filterwarnings('ignore')
    import subprocess
    from subprocess import DEVNULL, STDOUT, check_call  
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)

# Aim: Find activate prophages using Propagate
# Note: This script should be run within conda env "Propagate" by "conda activate Propagate"


# Step 1 Store IMGID set
IMGID = set()
with open('TYMEFLIES_metagenome_info.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('IMG'):
            tmp = line.split('\t')
            IMGID.add(tmp[0])


# Step 2 Run Propagate for individual metagenomes
propagate_cmd = []
for IMG in IMGID:
    # Test if the inputs are present
    if os.path.isfile(f"{IMG}/VIBRANT_{IMG}.a.v2.min2000/prophage_coordinates.txt"):
        if os.path.isfile(f"{IMG}/{IMG}.id97.bam") and os.path.isfile(f"{IMG}/{IMG}.a.fna"):
            each_propagate_cmd = f"Propagate -f {IMG}/{IMG}.a.fna -b {IMG}/{IMG}.id97.bam -v {IMG}/VIBRANT_{IMG}.a.v2.min2000/prophage_coordinates.txt -o {IMG}/PropagAtE_result --clean -t 1"
            propagate_cmd.append(each_propagate_cmd)
           
n = 10 # The number of parallel processes
for j in range(max(int(len(propagate_cmd)/n + 1), 1)):
    procs = [subprocess.Popen(i, shell=True, stdout=DEVNULL) for i in propagate_cmd[j*n: min((j+1)*n, len(propagate_cmd))] ]
    for p in procs:
        p.wait()              

           