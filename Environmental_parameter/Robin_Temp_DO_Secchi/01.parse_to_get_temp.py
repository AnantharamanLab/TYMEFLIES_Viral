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
    
# Aim: Parse to get temperature results


# Step 1 Store the input table 
input_table = 'temp.tsv'
date2temp_list = defaultdict(list)
with open(input_table, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('Sample'):
            items = line.split('\t')
            if int(items[1]) >= 2000 and  int(items[1]) <= 2019:
                date = items[1] + '-' + str(items[2]).zfill(2) + '-' + str(items[3]).zfill(2)
                temp = str(items[34])
                if temp != 'NA':
                    date2temp_list[date].append(temp) 
                    
                    
# Step 2 Get the temp for temp list
date2temp = {} # date => temp
for date in sorted(date2temp_list.keys()):
    temp_list = date2temp_list[date]
    # Convert strings to floats
    temp_list = [float(temp) for temp in temp_list]
    # Calculate the mean
    mean_temp = sum(temp_list) / len(temp_list)
    date2temp[date] = str(mean_temp)
    
    
# Step 3 Write down the result
f = open('temp.parsed.tsv', 'w') 
for date in sorted(date2temp.keys()):  
    line = date + '\t' + date2temp[date]
    f.write(line + '\n')
    
        
            
            
            



