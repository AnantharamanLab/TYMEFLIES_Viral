import matplotlib.pyplot as plt
import matplotlib.pyplot as pltt
import numpy as np
from collections import defaultdict
from pathlib import Path
from glob import glob
import statistics

def remove_nan_and_calculate_statistics(data_list):
    # Convert 'nan' strings to actual NaN values (np.nan) and remove them from the list
    valid_data = [float(value) for value in data_list if str(value) != 'nan']
    
    # Calculate mean and std of valid data (at least three values required)
    if len(valid_data) >= 3:
        mean = np.mean(valid_data)
        std = np.std(valid_data)
        return mean, std
    else:
        return None, None

# Pre-defined parameters
x_mean = np.arange(-65, 150, 5) # Define mean line x points


# Step 1 Store interpolation results for the first figure of Methylocystis_MAG
## Store the highest y1 value mean
highest_y1_values = []
with open("Methylocystis_virus_n_MAG/Methylocystis_MAG.year2highest_y_value.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y1_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y1_value != 'nan':
            highest_y1_values.append(float(highest_y1_value))

## Calculate the mean of highest y1 values
highest_y1_value_mean = statistics.mean(highest_y1_values)

## Store interpolation results
interpolation_result_files1 = glob("Methylocystis_virus_n_MAG/*.Methylocystis_MAG.interpolation_results.txt")
x_data1 = x_mean
year2x2y1 = defaultdict(dict) # year => x => y
for interpolation_result_file in interpolation_result_files1:
    year = Path(interpolation_result_file).stem.split('_')[0]
    y_list = []
    with open(interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            x, y = line.split('\t')[0], line.split('\t')[1]
            y_list.append(y)
    year2x2y1[year]['y'] = y_list

## Prepare data for plotting for the first figure
y_mean_values1 = []
y_std_values1 = []

## Initialize a dictionary to store all y values for each x point
x2y_values1 = defaultdict(list)

### Process each year's data for the first figure
for year, data in year2x2y1.items():
    y_data = data['y']  # Get the y values for the current year

    # Process each x point's data
    for x_val, y_val in zip(x_data1, y_data):
        x2y_values1[x_val].append(y_val)

### Calculate statistics for each x point
for x_val, y_data in x2y_values1.items():
    y_mean, y_std = remove_nan_and_calculate_statistics(y_data)

    # Append mean values and std values to the list
    y_mean_values1.append(y_mean * highest_y1_value_mean / 100 if y_mean is not None else np.nan)
    y_std_values1.append(y_std * highest_y1_value_mean / 100 if y_std is not None else np.nan)
    

# Step 2 Store interpolation results for the second figure of pmoC_containing_Methylocystis_viral_gn
## Store the highest y2 value mean
highest_y2_values = []
with open("Methylocystis_virus_n_MAG/pmoC_containing_Methylocystis_viral_gn.year2highest_y_value.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y2_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y2_value != 'nan':
            highest_y2_values.append(float(highest_y2_value))

## Calculate the mean of highest y2 values
highest_y2_value_mean = statistics.mean(highest_y2_values)

## Store interpolation results
interpolation_result_files2 = glob("Methylocystis_virus_n_MAG/*.pmoC_containing_Methylocystis_viral_gn.interpolation_results.txt")
x_data2 = x_mean
year2x2y2 = defaultdict(dict) # year => x => y
for interpolation_result_file in interpolation_result_files2:
    year = Path(interpolation_result_file).stem.split('_')[0]
    y_list = []
    with open(interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            x, y = line.split('\t')[0], line.split('\t')[1]
            y_list.append(y)
    year2x2y2[year]['y'] = y_list

## Prepare data for plotting for the second figure
y_mean_values2 = []
y_std_values2 = []

## Initialize a dictionary to store all y values for each x point
x2y_values2 = defaultdict(list)

### Process each year's data for the second figure
for year, data in year2x2y2.items():
    y_data = data['y']  # Get the y values for the current year

    # Process each x point's data
    for x_val, y_val in zip(x_data2, y_data):
        x2y_values2[x_val].append(y_val)

### Calculate statistics for each x point
for x_val, y_data in x2y_values2.items():
    y_mean, y_std = remove_nan_and_calculate_statistics(y_data)

    # Append mean values and std values to the list
    y_mean_values2.append(y_mean * highest_y2_value_mean / 100 if y_mean is not None else np.nan)
    y_std_values2.append(y_std * highest_y2_value_mean / 100 if y_std is not None else np.nan)

    
# Step 3 Store interpolation results for the third figure of UBA6136_MAG
## Store the highest y3 value mean
highest_y3_values = []
with open("UBA6136_virus_n_MAG/UBA6136_MAG.year2highest_y_value.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y3_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y3_value != 'nan':
            highest_y3_values.append(float(highest_y3_value))

## Calculate the mean of highest y3 values
highest_y3_value_mean = statistics.mean(highest_y3_values)

## Store interpolation results
interpolation_result_files3 = glob("UBA6136_virus_n_MAG/*.UBA6136_MAG.interpolation_results.txt")
x_data3 = x_mean
year2x2y3 = defaultdict(dict) # year => x => y
for interpolation_result_file in interpolation_result_files3:
    year = Path(interpolation_result_file).stem.split('_')[0]
    y_list = []
    with open(interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            x, y = line.split('\t')[0], line.split('\t')[1]
            y_list.append(y)
    year2x2y3[year]['y'] = y_list

## Prepare data for plotting for the third figure
y_mean_values3 = []
y_std_values3 = []

## Initialize a dictionary to store all y values for each x point
x2y_values3 = defaultdict(list)

### Process each year's data for the third figure
for year, data in year2x2y3.items():
    y_data = data['y']  # Get the y values for the current year

    # Process each x point's data
    for x_val, y_val in zip(x_data3, y_data):
        x2y_values3[x_val].append(y_val)

### Calculate statistics for each x point
for x_val, y_data in x2y_values3.items():
    y_mean, y_std = remove_nan_and_calculate_statistics(y_data)

    # Append mean values and std values to the list
    y_mean_values3.append(y_mean * highest_y3_value_mean / 100 if y_mean is not None else np.nan)
    y_std_values3.append(y_std * highest_y3_value_mean / 100 if y_std is not None else np.nan)
    

# Step 4 Store interpolation results for the fourth figure of pmoC_containing_UBA6136_viral_gn
## Store the highest y4 value mean
highest_y4_values = []
with open("UBA6136_virus_n_MAG/pmoC_containing_UBA6136_viral_gn.year2highest_y_value.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y4_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y4_value != 'nan':
            highest_y4_values.append(float(highest_y4_value))

## Calculate the mean of highest y4 values
highest_y4_value_mean = statistics.mean(highest_y4_values)

## Store interpolation results
interpolation_result_files4 = glob("UBA6136_virus_n_MAG/*.pmoC_containing_UBA6136_viral_gn.interpolation_results.txt")
x_data4 = x_mean
year2x2y4 = defaultdict(dict) # year => x => y
for interpolation_result_file in interpolation_result_files4:
    year = Path(interpolation_result_file).stem.split('_')[0]
    y_list = []
    with open(interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            x, y = line.split('\t')[0], line.split('\t')[1]
            y_list.append(y)
    year2x2y4[year]['y'] = y_list

## Prepare data for plotting for the fourth figure
y_mean_values4 = []
y_std_values4 = []

## Initialize a dictionary to store all y values for each x point
x2y_values4 = defaultdict(list)

### Process each year's data for the fourth figure
for year, data in year2x2y4.items():
    y_data = data['y']  # Get the y values for the current year

    # Process each x point's data
    for x_val, y_val in zip(x_data4, y_data):
        x2y_values4[x_val].append(y_val)

### Calculate statistics for each x point
for x_val, y_data in x2y_values4.items():
    y_mean, y_std = remove_nan_and_calculate_statistics(y_data)

    # Append mean values and std values to the list
    y_mean_values4.append(y_mean * highest_y4_value_mean / 100 if y_mean is not None else np.nan)
    y_std_values4.append(y_std * highest_y4_value_mean / 100 if y_std is not None else np.nan)


# Step 5 Store interpolation results for the fifth figure of Methylomonas_MAG
## Store the highest y5 value mean
highest_y5_values = []
with open("Methylomonas_virus_n_MAG/Methylomonas_MAG.year2highest_y_value.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y5_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y5_value != 'nan':
            highest_y5_values.append(float(highest_y5_value))

## Calculate the mean of highest y5 values
highest_y5_value_mean = statistics.mean(highest_y5_values)

## Store interpolation results
interpolation_result_files5 = glob("Methylomonas_virus_n_MAG/*.Methylomonas_MAG.interpolation_results.txt")
x_data5 = x_mean
year2x2y5 = defaultdict(dict) # year => x => y
for interpolation_result_file in interpolation_result_files5:
    year = Path(interpolation_result_file).stem.split('_')[0]
    y_list = []
    with open(interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            x, y = line.split('\t')[0], line.split('\t')[1]
            y_list.append(y)
    year2x2y5[year]['y'] = y_list

## Prepare data for plotting for the fifth figure
y_mean_values5 = []
y_std_values5 = []

## Initialize a dictionary to store all y values for each x point
x2y_values5 = defaultdict(list)

### Process each year's data for the fifth figure
for year, data in year2x2y5.items():
    y_data = data['y']  # Get the y values for the current year

    # Process each x point's data
    for x_val, y_val in zip(x_data5, y_data):
        x2y_values5[x_val].append(y_val)

### Calculate statistics for each x point
for x_val, y_data in x2y_values5.items():
    y_mean, y_std = remove_nan_and_calculate_statistics(y_data)

    # Append mean values and std values to the list
    y_mean_values5.append(y_mean * highest_y5_value_mean / 100 if y_mean is not None else np.nan)
    y_std_values5.append(y_std * highest_y5_value_mean / 100 if y_std is not None else np.nan)
    

# Step 6 Store interpolation results for the sixth figure of pmoC_containing_Methylomonas_viral_gn
## Store the highest y6 value mean
highest_y6_values = []
with open("Methylomonas_virus_n_MAG/pmoC_containing_Methylomonas_viral_gn.year2highest_y_value.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y6_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y6_value != 'nan':
            highest_y6_values.append(float(highest_y6_value))

## Calculate the mean of highest y6 values
highest_y6_value_mean = statistics.mean(highest_y6_values)

## Store interpolation results
interpolation_result_files6 = glob("Methylomonas_virus_n_MAG/*.pmoC_containing_Methylomonas_viral_gn.interpolation_results.txt")
x_data6 = x_mean
year2x2y6 = defaultdict(dict) # year => x => y
for interpolation_result_file in interpolation_result_files6:
    year = Path(interpolation_result_file).stem.split('_')[0]
    y_list = []
    with open(interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            x, y = line.split('\t')[0], line.split('\t')[1]
            y_list.append(y)
    year2x2y6[year]['y'] = y_list

## Prepare data for plotting for the sixth figure
y_mean_values6 = []
y_std_values6 = []

## Initialize a dictionary to store all y values for each x point
x2y_values6 = defaultdict(list)

### Process each year's data for the sixth figure
for year, data in year2x2y6.items():
    y_data = data['y']  # Get the y values for the current year

    # Process each x point's data
    for x_val, y_val in zip(x_data6, y_data):
        x2y_values6[x_val].append(y_val)

### Calculate statistics for each x point
for x_val, y_data in x2y_values6.items():
    y_mean, y_std = remove_nan_and_calculate_statistics(y_data)

    # Append mean values and std values to the list
    y_mean_values6.append(y_mean * highest_y6_value_mean / 100 if y_mean is not None else np.nan)
    y_std_values6.append(y_std * highest_y6_value_mean / 100 if y_std is not None else np.nan)


# Step 7 Store interpolation results for the seventh figure of UBA10906_MAG
## Store the highest y7 value mean
highest_y7_values = []
with open("UBA10906_virus_n_MAG/UBA10906_MAG.year2highest_y_value.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y7_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y7_value != 'nan':
            highest_y7_values.append(float(highest_y7_value))

## Calculate the mean of highest y7 values
highest_y7_value_mean = statistics.mean(highest_y7_values)

## Store interpolation results
interpolation_result_files7 = glob("UBA10906_virus_n_MAG/*.UBA10906_MAG.interpolation_results.txt")
x_data7 = x_mean
year2x2y7 = defaultdict(dict) # year => x => y
for interpolation_result_file in interpolation_result_files7:
    year = Path(interpolation_result_file).stem.split('_')[0]
    y_list = []
    with open(interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            x, y = line.split('\t')[0], line.split('\t')[1]
            y_list.append(y)
    year2x2y7[year]['y'] = y_list

## Prepare data for plotting for the seventh figure
y_mean_values7 = []
y_std_values7 = []

## Initialize a dictionary to store all y values for each x point
x2y_values7 = defaultdict(list)

### Process each year's data for the seventh figure
for year, data in year2x2y7.items():
    y_data = data['y']  # Get the y values for the current year

    # Process each x point's data
    for x_val, y_val in zip(x_data7, y_data):
        x2y_values7[x_val].append(y_val)

### Calculate statistics for each x point
for x_val, y_data in x2y_values7.items():
    y_mean, y_std = remove_nan_and_calculate_statistics(y_data)

    # Append mean values and std values to the list
    y_mean_values7.append(y_mean * highest_y7_value_mean / 100 if y_mean is not None else np.nan)
    y_std_values7.append(y_std * highest_y7_value_mean / 100 if y_std is not None else np.nan)
    

# Step 8 Store interpolation results for the eighth figure of pmoC_containing_UBA10906_viral_gn
## Store the highest y8 value mean
highest_y8_values = []
with open("UBA10906_virus_n_MAG/pmoC_containing_UBA10906_viral_gn.year2highest_y_value.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y8_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y8_value != 'nan':
            highest_y8_values.append(float(highest_y8_value))

## Calculate the mean of highest y8 values
highest_y8_value_mean = statistics.mean(highest_y8_values)

## Store interpolation results
interpolation_result_files8 = glob("UBA10906_virus_n_MAG/*.pmoC_containing_UBA10906_viral_gn.interpolation_results.txt")
x_data8 = x_mean
year2x2y8 = defaultdict(dict) # year => x => y
for interpolation_result_file in interpolation_result_files8:
    year = Path(interpolation_result_file).stem.split('_')[0]
    y_list = []
    with open(interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            x, y = line.split('\t')[0], line.split('\t')[1]
            y_list.append(y)
    year2x2y8[year]['y'] = y_list

## Prepare data for plotting for the eighth figure
y_mean_values8 = []
y_std_values8 = []

## Initialize a dictionary to store all y values for each x point
x2y_values8 = defaultdict(list)

### Process each year's data for the eighth figure
for year, data in year2x2y8.items():
    y_data = data['y']  # Get the y values for the current year

    # Process each x point's data
    for x_val, y_val in zip(x_data8, y_data):
        x2y_values8[x_val].append(y_val)

### Calculate statistics for each x point
for x_val, y_data in x2y_values8.items():
    y_mean, y_std = remove_nan_and_calculate_statistics(y_data)

    # Append mean values and std values to the list
    y_mean_values8.append(y_mean * highest_y8_value_mean / 100 if y_mean is not None else np.nan)
    y_std_values8.append(y_std * highest_y8_value_mean / 100 if y_std is not None else np.nan)
    
    
# Step 9 Store interpolation results for the nineth figure of no_pmoC_containing_Methylocystis_viral_gn
## Store the highest y9 value mean
highest_y9_values = []
with open("Methylocystis_virus_n_MAG/no_pmoC_containing_Methylocystis_viral_gn.year2highest_y_value.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y9_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y9_value != 'nan':
            highest_y9_values.append(float(highest_y9_value))

## Calculate the mean of highest y9 values
highest_y9_value_mean = statistics.mean(highest_y9_values)

## Store interpolation results
interpolation_result_files9 = glob("Methylocystis_virus_n_MAG/*.no_pmoC_containing_Methylocystis_viral_gn.interpolation_results.txt")
x_data9 = x_mean
year2x2y9 = defaultdict(dict) # year => x => y
for interpolation_result_file in interpolation_result_files9:
    year = Path(interpolation_result_file).stem.split('_')[0]
    y_list = []
    with open(interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            x, y = line.split('\t')[0], line.split('\t')[1]
            y_list.append(y)
    year2x2y9[year]['y'] = y_list

## Prepare data for plotting for the nineth figure
y_mean_values9 = []
y_std_values9 = []

## Initialize a dictionary to store all y values for each x point
x2y_values9 = defaultdict(list)

### Process each year's data for the nineth figure
for year, data in year2x2y9.items():
    y_data = data['y']  # Get the y values for the current year

    # Process each x point's data
    for x_val, y_val in zip(x_data9, y_data):
        x2y_values9[x_val].append(y_val)

### Calculate statistics for each x point
for x_val, y_data in x2y_values9.items():
    y_mean, y_std = remove_nan_and_calculate_statistics(y_data)

    # Append mean values and std values to the list
    y_mean_values9.append(y_mean * highest_y9_value_mean / 100 if y_mean is not None else np.nan)
    y_std_values9.append(y_std * highest_y9_value_mean / 100 if y_std is not None else np.nan)  


# Step 10 Store interpolation results for the 10th figure of no_pmoC_containing_UBA6136_viral_gn
## Store the highest y10 value mean
highest_y10_values = []
with open("UBA6136_virus_n_MAG/no_pmoC_containing_UBA6136_viral_gn.year2highest_y_value.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y10_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y10_value != 'nan':
            highest_y10_values.append(float(highest_y10_value))

## Calculate the mean of highest y10 values
highest_y10_value_mean = statistics.mean(highest_y10_values)

## Store interpolation results
interpolation_result_files10 = glob("UBA6136_virus_n_MAG/*.no_pmoC_containing_UBA6136_viral_gn.interpolation_results.txt")
x_data10 = x_mean
year2x2y10 = defaultdict(dict) # year => x => y
for interpolation_result_file in interpolation_result_files10:
    year = Path(interpolation_result_file).stem.split('_')[0]
    y_list = []
    with open(interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            x, y = line.split('\t')[0], line.split('\t')[1]
            y_list.append(y)
    year2x2y10[year]['y'] = y_list

## Prepare data for plotting for the 10th figure
y_mean_values10 = []
y_std_values10 = []

## Initialize a dictionary to store all y values for each x point
x2y_values10 = defaultdict(list)

### Process each year's data for the 10th figure
for year, data in year2x2y10.items():
    y_data = data['y']  # Get the y values for the current year

    # Process each x point's data
    for x_val, y_val in zip(x_data10, y_data):
        x2y_values10[x_val].append(y_val)

### Calculate statistics for each x point
for x_val, y_data in x2y_values10.items():
    y_mean, y_std = remove_nan_and_calculate_statistics(y_data)

    # Append mean values and std values to the list
    y_mean_values10.append(y_mean * highest_y10_value_mean / 100 if y_mean is not None else np.nan)
    y_std_values10.append(y_std * highest_y10_value_mean / 100 if y_std is not None else np.nan)  


# Step 11 Store interpolation results for the 11th figure of no_pmoC_containing_Methylomonas_viral_gn
## Store the highest y11 value mean
highest_y11_values = []
with open("Methylomonas_virus_n_MAG/no_pmoC_containing_Methylomonas_viral_gn.year2highest_y_value.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y11_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y11_value != 'nan':
            highest_y11_values.append(float(highest_y11_value))

## Calculate the mean of highest y11 values
highest_y11_value_mean = statistics.mean(highest_y11_values)

## Store interpolation results
interpolation_result_files11 = glob("Methylomonas_virus_n_MAG/*.no_pmoC_containing_Methylomonas_viral_gn.interpolation_results.txt")
x_data11 = x_mean
year2x2y11 = defaultdict(dict) # year => x => y
for interpolation_result_file in interpolation_result_files11:
    year = Path(interpolation_result_file).stem.split('_')[0]
    y_list = []
    with open(interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            x, y = line.split('\t')[0], line.split('\t')[1]
            y_list.append(y)
    year2x2y11[year]['y'] = y_list

## Prepare data for plotting for the 11th figure
y_mean_values11 = []
y_std_values11 = []

## Initialize a dictionary to store all y values for each x point
x2y_values11 = defaultdict(list)

### Process each year's data for the 11th figure
for year, data in year2x2y11.items():
    y_data = data['y']  # Get the y values for the current year

    # Process each x point's data
    for x_val, y_val in zip(x_data11, y_data):
        x2y_values11[x_val].append(y_val)

### Calculate statistics for each x point
for x_val, y_data in x2y_values11.items():
    y_mean, y_std = remove_nan_and_calculate_statistics(y_data)

    # Append mean values and std values to the list
    y_mean_values11.append(y_mean * highest_y11_value_mean / 100 if y_mean is not None else np.nan)
    y_std_values11.append(y_std * highest_y11_value_mean / 100 if y_std is not None else np.nan)    
    
    
# Step 12 Store interpolation results for the 12th figure of no_pmoC_containing_UBA10906_viral_gn
## Store the highest y12 value mean
highest_y12_values = []
with open("UBA10906_virus_n_MAG/no_pmoC_containing_UBA10906_viral_gn.year2highest_y_value.txt", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y12_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y12_value != 'nan':
            highest_y12_values.append(float(highest_y12_value))

## Calculate the mean of highest y12 values
highest_y12_value_mean = statistics.mean(highest_y12_values)

## Store interpolation results
interpolation_result_files12 = glob("UBA10906_virus_n_MAG/*.no_pmoC_containing_UBA10906_viral_gn.interpolation_results.txt")
x_data12 = x_mean
year2x2y12 = defaultdict(dict) # year => x => y
for interpolation_result_file in interpolation_result_files12:
    year = Path(interpolation_result_file).stem.split('_')[0]
    y_list = []
    with open(interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            x, y = line.split('\t')[0], line.split('\t')[1]
            y_list.append(y)
    year2x2y12[year]['y'] = y_list

## Prepare data for plotting for the 12th figure
y_mean_values12 = []
y_std_values12 = []

## Initialize a dictionary to store all y values for each x point
x2y_values12 = defaultdict(list)

### Process each year's data for the 12th figure
for year, data in year2x2y12.items():
    y_data = data['y']  # Get the y values for the current year

    # Process each x point's data
    for x_val, y_val in zip(x_data12, y_data):
        x2y_values12[x_val].append(y_val)

### Calculate statistics for each x point
for x_val, y_data in x2y_values12.items():
    y_mean, y_std = remove_nan_and_calculate_statistics(y_data)

    # Append mean values and std values to the list
    y_mean_values12.append(y_mean * highest_y12_value_mean / 100 if y_mean is not None else np.nan)
    y_std_values12.append(y_std * highest_y12_value_mean / 100 if y_std is not None else np.nan)      


# Step 13 Plotting the eight mean lines together in one figure
plt.figure(figsize=(5, 2.5))
plt.xlim(-65, 145)
plt.ylim(0, 60)  # Adjust the ylim to accommodate the first Axis range

## Plot the Methylocystis_MAG
plt.plot(x_data1, y_mean_values1, color='#0099CC', linestyle='solid', linewidth=0.5)
## Plot the pmoC_containing_Methylocystis_viral_gn
plt.plot(x_data2, y_mean_values2, color='#0099CC', linestyle='dotted', linewidth=0.5)

## Plot the UBA6136_MAG
plt.plot(x_data3, y_mean_values3, color='#336666', linestyle='solid', linewidth=0.5)
## Plot the pmoC_containing_UBA6136_viral_gn
plt.plot(x_data4, y_mean_values4, color='#336666', linestyle='dotted', linewidth=0.5)

## Plot the Methylomonas_MAG
plt.plot(x_data5, y_mean_values5, color='#FF6666', linestyle='solid', linewidth=0.5)
## Plot the pmoC_containing_Methylomonas_viral_gn
plt.plot(x_data6, y_mean_values6, color='#FF6666', linestyle='dotted', linewidth=0.5)

## Plot the UBA10906_MAG
plt.plot(x_data7, y_mean_values7, color='#F0E68C', linestyle='solid', linewidth=0.5)
## Plot the pmoC_containing_UBA10906_viral_gn
plt.plot(x_data8, y_mean_values8, color='#F0E68C', linestyle='dotted', linewidth=0.5)

## Create a twin Axes for the other four mean lines
ax2 = plt.gca().twinx()
ax2.set_ylim(0, 800)  # Set the y-axis limits for the fifth mean line

## Plot the no_pmoC_containing_Methylocystis_viral_gn
ax2.plot(x_data9, y_mean_values9, color = '#0099CC', linestyle='dashed', linewidth=0.5)
## Plot the no_pmoC_containing_UBA6136_viral_gn
ax2.plot(x_data10, y_mean_values10, color = '#336666', linestyle='dashed', linewidth=0.5)
## Plot the no_pmoC_containing_Methylomonas_viral_gn
ax2.plot(x_data11, y_mean_values11, color = '#FF6666', linestyle='dashed', linewidth=0.5)
## Plot the no_pmoC_containing_UBA10906_viral_gn
ax2.plot(x_data12, y_mean_values12, color = '#F0E68C', linestyle='dashed', linewidth=0.5)

## Combine the legends for both twin Axes mean lines
lines, labels = plt.gca().get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
plt.legend(lines + lines2, labels + labels2)

plt.grid(True)

## Save the combined figure
plt.savefig('Methanotroph.combined.mean_lines.svg', format = 'svg')
