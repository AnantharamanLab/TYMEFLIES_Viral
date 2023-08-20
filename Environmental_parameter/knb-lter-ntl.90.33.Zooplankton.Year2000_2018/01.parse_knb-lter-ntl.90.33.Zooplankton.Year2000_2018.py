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
    
# Aim: Parse to get all Zooplankton results grouped by date


def change_date_format(date_str):
    # Convert the input string to a datetime object
    date_obj = datetime.datetime.strptime(date_str, '%Y-%m-%d')    
    # Format the datetime object in the desired format
    new_date_str = date_obj.strftime('%Y/%-m/%d')
    return new_date_str
    
    
# Step 1 Store the table 
input_table = 'TYMEFLIES_viral.EnvParameter.knb-lter-ntl.90.33.Zooplankton.Year2000_2018.txt'
line_list = [] # Store each line as a list
with open(input_table, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('Date'):
            line_list.append(line)
            

# Step 2 Get the date2density dict
date2density = defaultdict(int)  # date => density
for line in line_list:
    items = line.split('\t')
    date, density = items[1], int(items[3])    
    if date not in date2density:
        date2density[date] = density
    else:
        date2density[date] += density
    
    
# Step 3 Write down the results
output_file = 'TYMEFLIES_viral.EnvParameter.knb-lter-ntl.90.33.Zooplankton.Year2000_2018.parsed.txt'
with open(output_file, 'w') as file:
    for date in date2density:
        date1 = change_date_format(date)
        density = date2density[date]
        line = date1 + '\t' + date + '\t' + str(density)  
        file.write(line + '\n')        