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
output_directory = "UBA10906_virus_n_MAG"  # Replace with your desired output directory
x_mean = np.arange(-65, 150, 5) # Define mean line x points


# Step 1 Store interpolation results for the first figure of UBA10906_MAG
## Store the highest y1 value mean
highest_y1_values = []
with open(f"{Path(output_directory) / 'UBA10906_MAG.year2highest_y_value.txt'}", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y1_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y1_value != 'nan':
            highest_y1_values.append(float(highest_y1_value))

## Calculate the mean of highest y1 values
highest_y1_value_mean = statistics.mean(highest_y1_values)

## Store interpolation results
interpolation_result_files1 = glob(f"{Path(output_directory) / '*.UBA10906_MAG.interpolation_results.txt'}")
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
    

# Step 2 Store interpolation results for the second figure of pmoC_containing_UBA10906_viral_gn
## Store the highest y2 value mean
highest_y2_values = []
with open(f"{Path(output_directory) / 'pmoC_containing_UBA10906_viral_gn.year2highest_y_value.txt'}", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y2_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y2_value != 'nan':
            highest_y2_values.append(float(highest_y2_value))

## Calculate the mean of highest y2 values
highest_y2_value_mean = statistics.mean(highest_y2_values)

## Store interpolation results
interpolation_result_files2 = glob(f"{Path(output_directory) / '*.pmoC_containing_UBA10906_viral_gn.interpolation_results.txt'}")
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

    
# Step 3 Store interpolation results for the third figure of no_pmoC_containing_UBA10906_viral_gn
## Store the highest y3 value mean
highest_y3_values = []
with open(f"{Path(output_directory) / 'no_pmoC_containing_UBA10906_viral_gn.year2highest_y_value.txt'}", 'r') as lines:
    for line in lines:
        line = line.rstrip('\n')
        year, highest_y3_value = line.split('\t')[0], line.split('\t')[1]
        if highest_y3_value != 'nan':
            highest_y3_values.append(float(highest_y3_value))

## Calculate the mean of highest y3 values
highest_y3_value_mean = statistics.mean(highest_y3_values)

## Store interpolation results
interpolation_result_files3 = glob(f"{Path(output_directory) / '*.no_pmoC_containing_UBA10906_viral_gn.interpolation_results.txt'}")
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


# Step 4 Plotting the two mean lines together in one figure
plt.figure(figsize=(5, 2.5))
plt.xlim(-65, 145)
plt.ylim(0, 50)  # Adjust the ylim to accommodate the first line

## Plot the second and third mean line (pmoC_containing_UBA10906_viral_gn + no_pmoC_containing_UBA10906_viral_gn)
plt.plot(x_data2, y_mean_values2, color='#0099CC')
plt.fill_between(x_data2, np.subtract(y_mean_values2, y_std_values2), np.add(y_mean_values2, y_std_values2), alpha=0.2, color='#0099CC')

## Create a twin Axes for the third mean line
ax2 = plt.gca().twinx()
ax2.set_ylim(0, 800)  # Set the y-axis limits for the third mean line

## Plot the third mean line on the twin Axes
ax2.plot(x_data3, y_mean_values3, color = '#336666')
ax2.fill_between(x_data3, np.subtract(y_mean_values3, y_std_values3), np.add(y_mean_values3, y_std_values3), alpha=0.2, color = '#336666')

plt.xlabel('X Values')
plt.title('Mean and Std of Y values for each Year')

## Combine the legends for both mean lines
lines, labels = plt.gca().get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
plt.legend(lines + lines2, labels + labels2)

plt.grid(True)

## Save the combined figure
plt.savefig('UBA10906_viral_gn.combined.mean_lines.pdf')


# Step 5 Plotting the first mean line together in one figure
pltt.figure(figsize=(5, 2.5))
pltt.xlim(-65, 145)
pltt.ylim(0, 75)  # Adjust the ylim to accommodate the first line

## Plot the first mean line
pltt.plot(x_data1, y_mean_values1, color ='#FF6666')
pltt.fill_between(x_data1, np.subtract(y_mean_values1, y_std_values1), np.add(y_mean_values1, y_std_values1), alpha=0.2, color = '#FF6666')

pltt.xlabel('X Values')
pltt.title('Mean and Std of Y values for each Year')

pltt.grid(True)

## Save the combined figure
pltt.savefig('UBA10906_MAG.mean_line.v2.pdf')