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
    
# Aim: Parse all parameters by seasons


# Get the pre-set lists
## The 33 environmental parameter items
parameter_items = [
    "Temperature", "Dissolved oxygen", "Secchi depth", "Chlorophyll-a", 
    "Chloride", "Sulfate", "Calcium", "Magnesium", "Sodium", "Potassium", 
    "Iron", "Manganese", "pH", "Dissolved inorganic carbon", "Total inorganic carbon", 
    "Dissolved organic carbon", "Total organic carbon", "Nitrate plus nitrite", 
    "Ammonium", "Total phosphorus", "Soluble reactive phosphorus", 
    "Phytoplankton-Bacillariophyta", "Phytoplankton-Chlorophyta", "Phytoplankton-Chrysophyta", 
    "Phytoplankton-Cryptophyta", "Phytoplankton-Cyanophyta", "Phytoplankton-Euglenophyta", 
    "Phytoplankton-Haptophyta", "Phytoplankton-Miscellaneous", "Phytoplankton-Pyrrhophyta", 
    "Phytoplankton-Xanthophyta", "Phytoplankton-Total", "Zooplankton density"
]

## Get the year_season_list
season_list = ["Spring", "Clearwater", "Early Summer", "Late Summer", "Fall", "Ice-on"]
year_list = range(2000, 2020)
year_season_list = []
for year in year_list:
    for season in season_list:
        year_season = str(year) + '-' + season
        year_season_list.append(year_season)
        

def get_year_season(date_string):
    season = ''
    input_year, input_month, input_day = int(date_string.split('-')[0]), int(date_string.split('-')[1]), int(date_string.split('-')[2])
    input_date = datetime.date(input_year, input_month, input_day)
    start_date_Spring = datetime.date(int(Year2season2start_date[str(input_year)]['Spring'].split('-')[0]), int(Year2season2start_date[str(input_year)]['Spring'].split('-')[1]), int(Year2season2start_date[str(input_year)]['Spring'].split('-')[2]))
    start_date_Clearwater = datetime.date(int(Year2season2start_date[str(input_year)]['Clearwater'].split('-')[0]), int(Year2season2start_date[str(input_year)]['Clearwater'].split('-')[1]), int(Year2season2start_date[str(input_year)]['Clearwater'].split('-')[2]))
    start_date_Early_Summer = datetime.date(int(Year2season2start_date[str(input_year)]['Early Summer'].split('-')[0]), int(Year2season2start_date[str(input_year)]['Early Summer'].split('-')[1]), int(Year2season2start_date[str(input_year)]['Early Summer'].split('-')[2]))
    start_date_Late_Summer = datetime.date(int(Year2season2start_date[str(input_year)]['Late Summer'].split('-')[0]), int(Year2season2start_date[str(input_year)]['Late Summer'].split('-')[1]), int(Year2season2start_date[str(input_year)]['Late Summer'].split('-')[2]))
    start_date_Fall = datetime.date(int(Year2season2start_date[str(input_year)]['Fall'].split('-')[0]), int(Year2season2start_date[str(input_year)]['Fall'].split('-')[1]), int(Year2season2start_date[str(input_year)]['Fall'].split('-')[2]))
    start_date_Ice_on = datetime.date(int(Year2season2start_date[str(input_year)]['Ice-on'].split('-')[0]), int(Year2season2start_date[str(input_year)]['Ice-on'].split('-')[1]), int(Year2season2start_date[str(input_year)]['Ice-on'].split('-')[2]))
    
    start_date_year_start = datetime.date(input_year, 1, 1) 
    start_date_year_end = datetime.date(input_year, 12, 31) 
    
    if start_date_year_start <= input_date < start_date_Spring:
        season = 'Ice-on'
    elif start_date_Spring <= input_date < start_date_Clearwater:   
        season = 'Spring'
    elif start_date_Clearwater <= input_date < start_date_Early_Summer:
        season = 'Clearwater'
    elif start_date_Early_Summer <= input_date < start_date_Late_Summer:    
        season = 'Early Summer'
    elif start_date_Late_Summer <= input_date < start_date_Fall:    
        season = 'Late Summer'  
    elif start_date_Fall <= input_date < start_date_Ice_on:  
        season = 'Fall'  
    elif start_date_Ice_on <= input_date <= start_date_year_end:  
        season = 'Ice-on'  
    return str(input_year) + '-' + season
    

# Step 1 Store the "season_start_dates.txt"
Year2season2start_date = defaultdict(lambda: defaultdict(int)) # year => season => start_date
with open("/storage1/data11/TYMEFLIES_phage/season_start_dates.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('Year'):
            tmp = line.split('\t')
            year = str(tmp[0])
            Year2season2start_date[year]['Spring'] = tmp[1]
            Year2season2start_date[year]['Clearwater'] = tmp[2]
            Year2season2start_date[year]['Early Summer'] = tmp[3]
            Year2season2start_date[year]['Late Summer'] = tmp[4]
            Year2season2start_date[year]['Fall'] = tmp[5]
            Year2season2start_date[year]['Ice-on'] = tmp[6]
lines.close()


year_season2parameter2values = defaultdict(lambda: defaultdict(list))  # year_season => parameter => [values]
# Step 2 Parse temp
Header = [] # Store the header line
with open('Temp.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Date'):
            tmp = line.split('\t')
            Header = tmp
        else:
            tmp = line.split('\t')
            date = tmp[0]
            year_season = get_year_season(date)
            for i in range(1, len(tmp)):
                parameter = Header[i]
                value = tmp[i]
                year_season2parameter2values[year_season][parameter].append(value)


# Step 3 Parse DO
Header = [] # Store the header line
with open('DO.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Date'):
            tmp = line.split('\t')
            Header = tmp
        else:
            tmp = line.split('\t')
            date = tmp[0]
            year_season = get_year_season(date)
            for i in range(1, len(tmp)):
                parameter = Header[i]
                value = tmp[i]
                year_season2parameter2values[year_season][parameter].append(value)
                
              
# Step 4 Parse Secchi
Header = [] # Store the header line
with open('Secchi.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Date'):
            tmp = line.split('\t')
            Header = tmp
        else:
            tmp = line.split('\t')
            date = tmp[0]
            year_season = get_year_season(date)
            for i in range(1, len(tmp)):
                parameter = Header[i]
                value = tmp[i]
                year_season2parameter2values[year_season][parameter].append(value)
                
                
# Step 5 Parse Chl
Header = [] # Store the header line
with open('Chl.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Date'):
            tmp = line.split('\t')
            Header = tmp
        else:
            tmp = line.split('\t')
            date = tmp[0]
            year_season = get_year_season(date)
            for i in range(1, len(tmp)):
                parameter = Header[i]
                value = tmp[i]
                year_season2parameter2values[year_season][parameter].append(value)     


# Step 6 Parse Major ions
Header = [] # Store the header line
with open('Major_ions.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Date'):
            tmp = line.split('\t')
            Header = tmp
        else:
            tmp = line.split('\t')
            date = tmp[0]
            year_season = get_year_season(date)
            for i in range(1, len(tmp)):
                parameter = Header[i]
                value = tmp[i]
                year_season2parameter2values[year_season][parameter].append(value)  
                
                
# Step 7 Parse Chemical limnology
Header = [] # Store the header line
with open('Chemical_limnology.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Date'):
            tmp = line.split('\t')
            Header = tmp
        else:
            tmp = line.split('\t')
            date = tmp[0]
            year_season = get_year_season(date)
            for i in range(1, len(tmp)):
                parameter = Header[i]
                value = tmp[i]
                year_season2parameter2values[year_season][parameter].append(value)  


# Step 8 Parse Phytoplankton
Header = [] # Store the header line
with open('Phytoplankton.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Date'):
            tmp = line.split('\t')
            Header = tmp
        else:
            tmp = line.split('\t')
            date = tmp[0]
            year_season = get_year_season(date)
            for i in range(1, len(tmp)):
                parameter = Header[i]
                value = tmp[i]
                year_season2parameter2values[year_season][parameter].append(value)   


# Step 9 Parse Zooplankton
Header = [] # Store the header line
with open('Zooplankton.txt', 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('Date'):
            tmp = line.split('\t')
            Header = tmp
        else:
            tmp = line.split('\t')
            date = tmp[0]
            year_season = get_year_season(date)
            for i in range(1, len(tmp)):
                parameter = Header[i]
                value = tmp[i]
                year_season2parameter2values[year_season][parameter].append(value)                 
                
                                
# Step 10 Get the final summary result
## Step 10.1 Get year_season2parameter2value_mean dict
year_season2parameter2value_mean = defaultdict(lambda: defaultdict(dict))
for year_season in year_season_list:
    for parameter in parameter_items:
        value_mean = 'NA'
        if year_season in year_season2parameter2values and parameter in year_season2parameter2values[year_season]:
            values = year_season2parameter2values[year_season][parameter]
            values = [value for value in values if value != 'NA'] # Delete all 'NA's from the values list
            if len(values) > 0:
                # Convert strings to floats
                values = [float(value) for value in values]
                # Calculate the mean
                value_mean = sum(values) / len(values)
        year_season2parameter2value_mean[year_season][parameter] = value_mean
        
## Step 10.2 Write down the year_season2parameter2value_mean dict            
output_file = 'year_season2parameter2value_mean.txt'
with open(output_file, 'w') as file:
    ## Write table header
    header = '\t'.join(parameter_items)
    file.write('head\t' + header + '\n')
    
    ## Write table rows
    for year_season in year_season_list:
        row_values = [str(year_season2parameter2value_mean.get(year_season, {}).get(parameter, 'NA')) for parameter in parameter_items]
        row = year_season + '\t' + '\t'.join(row_values) + '\n'
        file.write(row)  
