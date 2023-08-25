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
    
# Aim: Parse iPHoP result to get the final host prediction result
# Host-genome predictions were ultimately determined using the following guideline: 
# the outcome (extracted from the result file: Host_prediction_to_genus_m90.csv) having the highest confidence score
#  was designated as the definitive outcome in cases where multiple prediction outcomes for a single virus are available. 


# Step 1 Store the result of Host_prediction_to_genus_m90.csv
virus2host_information = {} # virus => [virus, host_genus, confidence_score, list_of_methods]
with open("Host_prediction/iPHoP_result/Host_prediction_to_genus_m90.csv", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        if not line.startswith('Virus'):
            tmp = line.split(',')
            virus, host_genus, confidence_score, list_of_methods = tmp[0], tmp[2], tmp[3], tmp[4]
            if virus not in virus2host_information:
                virus2host_information[virus] = [virus, host_genus, confidence_score, list_of_methods]
            else: 
                if float(confidence_score) > float(virus2host_information[virus][2]):
                    virus2host_information[virus] = [virus, host_genus, confidence_score, list_of_methods]                 
   

# Step 2 Write down the parsed result
f = open('Host_prediction/iPHoP_result/iPHoP_result_parsed.txt', 'w')
for virus in virus2host_information:
    line = '\t'.join(virus2host_information[virus])
    f.write(line + '\n')
f.close()    
 

