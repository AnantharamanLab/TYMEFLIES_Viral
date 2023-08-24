#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    from collections import defaultdict
    import subprocess
    from subprocess import DEVNULL, STDOUT, check_call    
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: Calculate to get ahbD Early Summer AMG coverage


# Step 1 Store ahbD-containing viral gn coverage
viral_gn2earlysummer2cov = defaultdict(dict)
Header1 = []
EarlySummer_sets = set()
with open('/storage1/data11/TYMEFLIES_phage/MetaPop/AhbD_containing_viral_gn2year_season2cov.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Head'):
            tmp = line.split('\t')
            Header1 = tmp
        else:
            tmp = line.split('\t')
            viral_gn = tmp[0]
            for i in range(1, len(tmp)):
                if 'Early Summer' in Header1[i]:
                    earlysummer = Header1[i]
                    cov = tmp[i]    
                    viral_gn2earlysummer2cov[viral_gn][earlysummer] = cov
                    EarlySummer_sets.add(earlysummer)
                    

# Step 2 Store AMG cov ratio
AMG2earlysummer2cov_ratio = defaultdict(dict)
Header2 = []
with open('/storage1/data11/TYMEFLIES_phage/MetaPop/AhbD_AMG_gene2year_season2cov_ratio.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Head'):
            tmp = line.split('\t')
            Header2 = tmp   
        else:
            tmp = line.split('\t')  
            AMG = tmp[0]
            for i in range(1, len(tmp)):
                if 'Early Summer' in Header2[i]:
                    earlysummer = Header2[i]
                    cov_ratio = tmp[i]    
                    AMG2earlysummer2cov_ratio[AMG][earlysummer] = cov_ratio
            

# Step 3 Get ahbD AMG coverage and write down the result
AMG2earlysummer2cov = defaultdict(dict)
for AMG in AMG2earlysummer2cov_ratio:
    for earlysummer in EarlySummer_sets:
        cov_ratio = AMG2earlysummer2cov_ratio[AMG][earlysummer]
        viral_gn = AMG.split('__Ga', 1)[0]
        viral_gn_cov = viral_gn2earlysummer2cov[viral_gn][earlysummer]
        AMG_cov = float(viral_gn_cov) * float(cov_ratio)
        AMG2earlysummer2cov[AMG][earlysummer] = AMG_cov
        
f = open('/storage1/data11/TYMEFLIES_phage/MetaPop/AhbD_AMG_gene2earlysummer2cov.txt', 'w')
header = 'Head' + '\t' + '\t'.join(sorted(EarlySummer_sets)) + '\n'
f.write(header)
for AMG in sorted(AMG2earlysummer2cov.keys()):
    line = AMG + '\t'
    tmp = [] # Store all the cov for AMGs
    for earlysummer in sorted(EarlySummer_sets):
        tmp.append(str(AMG2earlysummer2cov[AMG][earlysummer]))    
    line += '\t'.join(tmp) + '\n'
    f.write(line)    
        



     
    