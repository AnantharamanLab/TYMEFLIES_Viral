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
    
# Aim: Parse to get Phytoplankton results grouped by divisions


def change_date_format(date_str):
    # Convert the input string to a datetime object
    date_obj = datetime.datetime.strptime(date_str, '%Y-%m-%d')    
    # Format the datetime object in the desired format
    new_date_str = date_obj.strftime('%Y/%-m/%d')
    return new_date_str
    
    
# Step 1 Store the table 
input_table = 'TYMEFLIES_viral.EnvParameter.knb-lter-ntl.88.31.Phytoplankton.Year2000_2019.txt'
line_list = [] # Store each line as a list
with open(input_table, 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('Date'):
            line_list.append(line)
            

# Step 2 Get the date2division2biomass_conc dict
date2division2biomass_conc = defaultdict(lambda: defaultdict(float))  # date => division => biomass_conc
division_set = set()

for line in line_list:
    items = line.split('\t')
    date, division, biomass_conc = items[1], items[2], float(items[4])
    division_set.add(division)
    
    if division not in date2division2biomass_conc[date]:
        date2division2biomass_conc[date][division] = biomass_conc
    else:
        date2division2biomass_conc[date][division] += biomass_conc
    
    
# Step 3 Write down the results
output_file = 'TYMEFLIES_viral.EnvParameter.knb-lter-ntl.88.31.Phytoplankton.Year2000_2019.parsed.txt'
with open(output_file, 'w') as file:
    ## Write table header
    header = '\t'.join(sorted(division_set))
    file.write('Date1\tDate\t' + header + '\n')
    
    ## Write table rows
    for date in sorted(date2division2biomass_conc.keys()):
        row_values = [str(date2division2biomass_conc.get(date, {}).get(division, '0')) for division in sorted(division_set)]
        date1 = change_date_format(date)
        row = date1 + '\t' + date + '\t' + '\t'.join(row_values) + '\n'
        file.write(row)     