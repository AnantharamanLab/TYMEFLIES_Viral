#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    import datetime
    from collections import defaultdict
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: Get date to season result and replace "TYMEFLIES_metagenome_info.txt"


def get_season(date_string):
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
    return season
    

# Step 1 Store the "season_start_dates.txt"
Year2season2start_date = defaultdict(lambda: defaultdict(int)) # year => season => start_date
with open("season_start_dates.txt", 'r') as lines:
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


# Step 2 Store the "TYMEFLIES_metagenome_info.txt" 
TYMEFLIES_metagenome_info = {} # IMG => [line_split_list]
header = ''
with open("TYMEFLIES_metagenome_info.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if line.startswith('IMG'):
            header = line
        else:
            line_split_list = line.split('\t')
            IMG, date_string, season = line_split_list[0], line_split_list[8], line_split_list[10] # An example for date_string is "YYYY-MM-DD"
            season = get_season(date_string)
            line_split_list[10] = season
            TYMEFLIES_metagenome_info[IMG] = line_split_list
lines.close()


# Step 3 Re-write the "TYMEFLIES_metagenome_info.txt"
f = open('TYMEFLIES_metagenome_info.new.txt', 'w')
f.write(header + '\n')
for IMG in TYMEFLIES_metagenome_info:
    line = '\t'.join(TYMEFLIES_metagenome_info[IMG]) + '\n'
    f.write(line)
f.close()    
            


            



                 