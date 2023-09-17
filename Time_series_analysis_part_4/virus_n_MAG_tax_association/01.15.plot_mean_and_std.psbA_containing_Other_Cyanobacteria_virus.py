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
        

# Step 1 Store the psbA-containing Other Cyanobacteria viral gn cov files for each year
psbA_containing_Other_Cyanobacteria_viral_gn_cov_files = glob("/storage1/data11/TYMEFLIES_phage/virus_n_MAG_tax_association/Other_Cyanobacteria_virus_n_MAG/*.psbA_containing_Other_Cyanobacteria_viral_gn_cov.txt")
year2psbA_containing_Other_Cyanobacteria_viral_gn_cov = defaultdict(dict)
for psbA_containing_Other_Cyanobacteria_viral_gn_cov_file in psbA_containing_Other_Cyanobacteria_viral_gn_cov_files:
    year = Path(psbA_containing_Other_Cyanobacteria_viral_gn_cov_file).stem.split('.')[0]    
    x_list = []
    y_list = []
    with open(psbA_containing_Other_Cyanobacteria_viral_gn_cov_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('earlysummer'):
                tmp = line.split('\t')
                x_list.append(float(tmp[0]))  # Convert x value to float
                y_list.append(float(tmp[1]))  # Convert y value to float
    year2psbA_containing_Other_Cyanobacteria_viral_gn_cov[year]['x'] = x_list 
    year2psbA_containing_Other_Cyanobacteria_viral_gn_cov[year]['y'] = y_list       


# Step 2 Generate the mean line
## Create interpolation functions for each year's dataset
interpolation_functions = {}
for year, data in year2psbA_containing_Other_Cyanobacteria_viral_gn_cov.items():
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
output_directory = "Other_Cyanobacteria_virus_n_MAG"  # Replace with your desired output directory

year2highest_y_value = {} # year => highest_y_value
for year, data in interpolation_results.items():
    output_filename = f"{year}.psbA_containing_Other_Cyanobacteria_viral_gn.interpolation_results.txt"
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

f = open('Other_Cyanobacteria_virus_n_MAG/psbA_containing_Other_Cyanobacteria_viral_gn.year2highest_y_value.txt', 'w')
for year in year2highest_y_value:
    f.write(f"{year}\t{year2highest_y_value[year]}\n")
f.close()    
    

# Step 5 Store interpolation results
interpolation_result_files = glob(f"{Path(output_directory) / '*.psbA_containing_Other_Cyanobacteria_viral_gn.interpolation_results.txt'}")
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
plt.savefig('psbA_containing_Other_Cyanobacteria_viral_gn.mean_line.pdf')


# Step 7 Get the late summer and fall dates bar plots
## Step 7.1 Store each year's Early Summer, Late Summer, Fall, Clearwater start datetime
year2earlysummer_start_day = {}
year2latesummer_start_day = {}
year2fall_start_day = {}
year2clearwater_start_day = {}
with open('/storage1/data11/TYMEFLIES_phage/season_start_dates.txt', 'r') as lines:
    for line in lines:
        if not line.startswith('Year'):
            line = line.rstrip('\n')
            year, earlysummer_start_day = line.split('\t')[0], line.split('\t')[3]
            latesummer_start_day, fall_start_day = line.split('\t')[4], line.split('\t')[5]
            clearwater_start_day = line.split('\t')[2]
            year2earlysummer_start_day[year] = earlysummer_start_day
            year2latesummer_start_day[year] = latesummer_start_day
            year2fall_start_day[year] = fall_start_day
            year2clearwater_start_day[year] = clearwater_start_day
            
## Step 7.2 Calculate the Late Summer, Fall, and Clearwater bars
latesummer_distance_list = []
fall_distance_list = []
clearwater_distance_list = []
for year in year2earlysummer_start_day:          
    earlysummer_start_day = year2earlysummer_start_day[year]
    latesummer_start_day = year2latesummer_start_day[year]
    fall_start_day = year2fall_start_day[year]
    clearwater_start_day = year2clearwater_start_day[year]
    date1 = datetime.strptime(earlysummer_start_day, "%Y-%m-%d")
    date2 = datetime.strptime(latesummer_start_day, "%Y-%m-%d")
    date3 = datetime.strptime(fall_start_day, "%Y-%m-%d")
    date4 = datetime.strptime(clearwater_start_day, "%Y-%m-%d")
    latesummer_distance_between_dates = date2 - date1
    fall_distance_between_dates = date3 - date1
    clearwater_distance_between_dates = date4 - date1
    latesummer_distance_list.append(latesummer_distance_between_dates.days) 
    fall_distance_list.append(fall_distance_between_dates.days) 
    clearwater_distance_list.append(clearwater_distance_between_dates.days) 

### Calculate average and standard deviation
latesummer_avg_distance = np.mean(latesummer_distance_list)
latesummer_std_distance = np.std(latesummer_distance_list)

### Create the plot
plt.figure(figsize=(8, 2))
plt.barh(0, latesummer_avg_distance, xerr=latesummer_std_distance, height=0.5, color='#6797bf', alpha=0.7, capsize=5)

### Set the y-axis limit
plt.ylim(-0.5, 0.5)

### Set the x-axis label
plt.xlabel('Distance')

### Set the plot title
plt.title('Average and Standard Deviation of Late Summer Distance')

### Save the plot
plt.savefig('Other_Cyanobacteria_virus_n_MAG.latesummer_distance_barplot.pdf')


### Calculate average and standard deviation
fall_avg_distance = np.mean(fall_distance_list)
fall_std_distance = np.std(fall_distance_list)

### Create the plot
plt.figure(figsize=(8, 2))
plt.barh(0, fall_avg_distance, xerr=fall_std_distance, height=0.5, color='#6797bf', alpha=0.7, capsize=5)

### Set the y-axis limit
plt.ylim(-0.5, 0.5)

### Set the x-axis label
plt.xlabel('Distance')

### Set the plot title
plt.title('Average and Standard Deviation of Fall Distance')

### Save the plot
plt.savefig('Other_Cyanobacteria_virus_n_MAG.fall_distance_barplot.pdf')


### Calculate average and standard deviation
clearwater_avg_distance = np.mean(clearwater_distance_list)
clearwater_std_distance = np.std(clearwater_distance_list)

### Create the plot
plt.figure(figsize=(8, 2))
plt.barh(0, clearwater_avg_distance, xerr=clearwater_std_distance, height=0.5, color='#6797bf', alpha=0.7, capsize=5)

### Set the y-axis limit
plt.ylim(-0.5, 0.5)

### Set the x-axis label
plt.xlabel('Distance')

### Set the plot title
plt.title('Average and Standard Deviation of Fall Distance')

### Save the plot
plt.savefig('Other_Cyanobacteria_virus_n_MAG.clearwater_distance_barplot.pdf')
