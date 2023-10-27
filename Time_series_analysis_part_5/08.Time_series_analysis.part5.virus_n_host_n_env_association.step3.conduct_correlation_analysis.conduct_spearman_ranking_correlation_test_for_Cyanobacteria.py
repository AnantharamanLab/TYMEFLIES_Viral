#!/usr/bin/env python

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings('ignore')
    import datetime
    from collections import defaultdict
    import statistics
    import subprocess
    import scipy.stats
    import statsmodels.stats.multitest as multi
    import numpy as np    
except Exception as e:
    sys.stderr.write(str(e) + '\n\n')
    exit(1)
    
# Aim: Conduct correlation analysis between each pair of virus/host and env parameters


def calculate_spearman_fisher_z(list1, list2, alpha=0.05):
    # Initialize empty lists to store filtered values
    filtered_list1 = []
    filtered_list2 = []

    # Iterate through both lists simultaneously and exclude positions where either list has "NA"
    for val1, val2 in zip(list1, list2):
        if val1 != "NA" and val2 != "NA":
            filtered_list1.append(val1)
            filtered_list2.append(val2)

    # Calculate Spearman's rank correlation coefficient
    correlation_coefficient, _ = scipy.stats.spearmanr(filtered_list1, filtered_list2)

    # Apply Fisher z-transformation
    fisher_z = 0.5 * (np.log(1 + correlation_coefficient) - np.log(1 - correlation_coefficient))

    # Calculate the standard error of the Fisher z-transformed correlation
    n = len(filtered_list1)
    fisher_z_se = 1 / np.sqrt(n - 3)

    # Calculate the z-statistic
    z_statistic = fisher_z / fisher_z_se

    # Calculate the p-value for the two-tailed test
    p_value = 2 * (1 - scipy.stats.norm.cdf(abs(z_statistic)))

    # Check if the p-value is significant at the specified alpha level
    is_significant = p_value < alpha

    return [correlation_coefficient, p_value, is_significant]



# Step 1 Store the environmental parameter and virus/host input table
row_name_list = [] # Store all the row names
row_name2list = {} # row_name => [list]

with open('virus_n_host_n_env_association/virus_n_host_n_env_association_input_table.txt', 'r') as infile:
    for line in infile:
        line = line.strip()
        if not line.startswith("Head"):
            tmp = line.split("\t")
            row_name = tmp[0]
            tmp_wo_row_name = tmp[1:]
            row_name2list[row_name] = tmp_wo_row_name
            row_name_list.append(row_name)
            

# Step 2 Get the pair of comparison and conduct the Spearman's rank correlation analysis
biology_entity_list = row_name_list[0:4]
env_list = row_name_list[4:]
biology_entity_list2env_list2result = defaultdict(dict)
for item1 in biology_entity_list:
    for item2 in env_list:
        list1 = row_name2list[item1]
        list2 = row_name2list[item2]
        correlation_coefficient, p_value, is_significant = calculate_spearman_fisher_z(list1, list2)
        biology_entity_list2env_list2result[item1][item2] = [correlation_coefficient, p_value, is_significant]
          
          
# Step 3: Write the results to an output file
output_filename = 'virus_n_host_n_env_association/spearman_correlation_results.txt'

with open(output_filename, 'w') as outfile:
    # Write the header row
    outfile.write('Biology Entity\tEnvironment\tCorrelation Coefficient\tP-Value\tis_significant\n')

    # Iterate through the dictionary and write the results
    for biology_entity, env_results in biology_entity_list2env_list2result.items():
        for environment, results in env_results.items():
            correlation_coefficient, p_value, is_significant = results
            # Format and write the row
            outfile.write(f'{biology_entity}\t{environment}\t{correlation_coefficient}\t{p_value}\t{is_significant}\n')

print(f'Results written to {output_filename}')