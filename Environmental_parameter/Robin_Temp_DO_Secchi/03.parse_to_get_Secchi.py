#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    warnings.filterwarnings("ignore")
    from pathlib import Path
    from glob import glob
    import datetime
    from collections import defaultdict
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)
    
# Aim: Parse to get Secchi results


# Step 1 Store the input table 
input_table = 'Secchi.tsv'
date2Secchi_list = defaultdict(list)
with open(input_table, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('Year'):
            items = line.split('\t')
            if int(items[0]) >= 2000 and  int(items[0]) <= 2019:
                date = items[0] + '-' + str(items[1]).zfill(2) + '-' + str(items[2]).zfill(2)
                Secchi = str(items[3])
                if Secchi != 'NA':
                    date2Secchi_list[date].append(Secchi) 
                    
                    
# Step 2 Get the Secchi for Secchi list
date2Secchi = {} # date => Secchi
for date in sorted(date2Secchi_list.keys()):
    Secchi_list = date2Secchi_list[date]
    # Convert strings to floats
    Secchi_list = [float(Secchi) for Secchi in Secchi_list]
    # Calculate the mean
    mean_Secchi = sum(Secchi_list) / len(Secchi_list)
    date2Secchi[date] = str(mean_Secchi)
    
    
# Step 3 Write down the result
f = open('Secchi.parsed.tsv', 'w') 
for date in sorted(date2Secchi.keys()):  
    line = date + '\t' + date2Secchi[date]
    f.write(line + '\n')
    
        
            
            
            



