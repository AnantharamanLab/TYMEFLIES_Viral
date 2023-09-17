import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from glob import glob
from pathlib import Path
from collections import defaultdict
from datetime import datetime  


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
        

# Step 1 Store the no_psbA_containing_Microcystis_viral_gn_cov files for each year
no_psbA_containing_Microcystis_viral_gn_cov_files = glob("/storage1/data11/TYMEFLIES_phage/virus_n_MAG_tax_association/Microcystis_virus_n_MAG/*.no_psbA_containing_Microcystis_viral_gn_cov.txt")
year2no_psbA_containing_Microcystis_viral_gn_cov = defaultdict(dict)
for no_psbA_containing_Microcystis_viral_gn_cov_file in no_psbA_containing_Microcystis_viral_gn_cov_files:
    year = Path(no_psbA_containing_Microcystis_viral_gn_cov_file).stem.split('.')[0]    
    x_list = []
    y_list = []
    with open(no_psbA_containing_Microcystis_viral_gn_cov_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('earlysummer'):
                tmp = line.split('\t')
                x_list.append(float(tmp[0]))  # Convert x value to float
                y_list.append(float(tmp[1]))  # Convert y value to float
    year2no_psbA_containing_Microcystis_viral_gn_cov[year]['x'] = x_list 
    year2no_psbA_containing_Microcystis_viral_gn_cov[year]['y'] = y_list       


# Step 2 Generate the mean line
## Create interpolation functions for each year's dataset
interpolation_functions = {}
for year, data in year2no_psbA_containing_Microcystis_viral_gn_cov.items():
    x = np.array(data['x'])  # Convert x values to numpy array
    y = np.array(data['y'])  # Convert y values to numpy array
    interpolation_functions[year] = interp1d(x, y, bounds_error=False)

## Define mean line x points
x_mean = np.arange(-45, 180, 5)


# Step 3 Evaluate interpolation functions at x_mean for each year
interpolation_results = defaultdict(dict)
for year, interp_func in interpolation_functions.items():
    y_mean = interp_func(x_mean)
    interpolation_results[year]['x_mean'] = x_mean
    interpolation_results[year]['y_mean'] = y_mean


# Step 4 Store the interpolation results in text files
output_directory = "Microcystis_virus_n_MAG"  # Replace with your desired output directory

year2highest_y_value = {} # year => highest_y_value
for year, data in interpolation_results.items():
    output_filename = f"{year}.no_psbA_containing_Microcystis_viral_gn.interpolation_results.txt"
    output_file_path = Path(output_directory) / output_filename

    x_mean_values = data['x_mean']
    y_mean_values = data['y_mean']

    # Find the highest value in y_mean_values
    highest_y_value = max(y_mean_values)
    year2highest_y_value[year] = highest_y_value
    
    with open(output_file_path, 'w') as output_file:
        for x_val, y_val in zip(x_mean_values, y_mean_values):
            # Calculate the percentage of y_val
            percentage_y_val = (y_val / highest_y_value) * 100

            # Print the modified y_val
            output_file.write(f"{x_val}\t{percentage_y_val:.2f}\n")

    print(f"Interpolation results for year {year} written to {output_file_path}")

f = open('Microcystis_virus_n_MAG/no_psbA_containing_Microcystis_viral_gn.year2highest_y_value.txt', 'w')
for year in year2highest_y_value:
    f.write(f"{year}\t{year2highest_y_value[year]}\n")
f.close()
    

# Step 5 Store interpolation results
interpolation_result_files = glob(f"{Path(output_directory) / '*.no_psbA_containing_Microcystis_viral_gn.interpolation_results.txt'}")
x_data = x_mean
year2x2y = defaultdict(dict) # year => x => y
for interpolation_result_file in interpolation_result_files:
    year = Path(interpolation_result_file).stem.split('_')[0]
    y_list = []
    with open(interpolation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            x, y = line.split('\t')[0], line.split('\t')[1]
            y_list.append(y)
    year2x2y[year]['y'] = y_list
     
     
# Step 6 Visualize the result
## Prepare data for plotting
years = []
y_mean_values = []
y_std_values = []

# Initialize a dictionary to store all y values for each x point
x2y_values = defaultdict(list)

## Process each year's data
for year, data in year2x2y.items():
    y_data = data['y']  # Get the y values for the current year

    # Process each x point's data
    for x_val, y_val in zip(x_data, y_data):
        x2y_values[x_val].append(y_val)

# Calculate statistics for each x point
for x_val, y_data in x2y_values.items():
    y_mean, y_std = remove_nan_and_calculate_statistics(y_data)

    # Append mean values and std values to the list
    y_mean_values.append(y_mean if y_mean is not None else np.nan)
    y_std_values.append(y_std if y_std is not None else np.nan) 


## Plotting the mean line and std line
plt.figure(figsize=(8, 5))
plt.xlim(-45, 180)
plt.ylim(0, 100)
plt.plot(x_data, y_mean_values, label='Mean of Y')
plt.fill_between(x_data, np.subtract(y_mean_values, y_std_values), np.add(y_mean_values, y_std_values), alpha=0.2)
plt.xlabel('X Values')
plt.ylabel('Y Values')
plt.title('Mean and Std of Y values for each Year')
plt.legend()
plt.grid(True)          
plt.savefig('no_psbA_containing_Microcystis_viral_gn.mean_line.pdf')