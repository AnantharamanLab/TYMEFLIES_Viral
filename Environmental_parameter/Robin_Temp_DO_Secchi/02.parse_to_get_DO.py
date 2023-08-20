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
    
# Aim: Parse to get DO results


# Step 1 Store the input table 
input_table = 'DO.tsv'
date2DO_list = defaultdict(list)
with open(input_table, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('Sample'):
            items = line.split('\t')
            if int(items[1]) >= 2000 and  int(items[1]) <= 2019:
                date = items[1] + '-' + str(items[2]).zfill(2) + '-' + str(items[3]).zfill(2)
                DO = str(items[36])
                if DO != 'NA':
                    date2DO_list[date].append(DO) 
                    
                    
# Step 2 Get the DO for DO list
date2DO = {} # date => DO
for date in sorted(date2DO_list.keys()):
    DO_list = date2DO_list[date]
    # Convert strings to floats
    DO_list = [float(DO) for DO in DO_list]
    # Calculate the mean
    mean_DO = sum(DO_list) / len(DO_list)
    date2DO[date] = str(mean_DO)
    
    
# Step 3 Write down the result
f = open('DO.parsed.tsv', 'w') 
for date in sorted(date2DO.keys()):  
    line = date + '\t' + date2DO[date]
    f.write(line + '\n')
    
        
            
            
            



