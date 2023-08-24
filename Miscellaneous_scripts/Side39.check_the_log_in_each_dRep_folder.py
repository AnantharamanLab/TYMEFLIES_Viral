#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    from collections import defaultdict
    from pathlib import Path
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
    
# Aim: Check the log file in each dRep folder


# Step 1 Store the dRep folder address dict and log file address dict
dRep_folder_addrs = glob("/storage1/data11/TYMEFLIES_phage/dRep_working_dir/Output.Genus_cluster*")
log_file_addrs = glob("/storage1/data11/TYMEFLIES_phage/dRep_working_dir/Output.Genus_cluster*/log/logger.log")
print(f"The number of dRep folders is {len(dRep_folder_addrs)}")
print(f"The number of log files is {len(log_file_addrs)}")

# Step 2 Check the log file
for log_file_addr in log_file_addrs:
    last_line_content = ''
    folder_name = log_file_addr.split('/')[5]
    search_string = 'Finished the dereplicate operation'
    with open(log_file_addr, 'r') as file:
        lines = file.readlines()
        last_line = lines[-1].strip()
        last_line_content = last_line
        
    if search_string in last_line_content:
        pass
    else:    
        print(f"{folder_name} was not fully finished")   