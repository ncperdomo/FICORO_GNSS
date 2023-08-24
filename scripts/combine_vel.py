import os
import sys
import pandas as pd
from itertools import combinations
import numpy as np
from math import sin, cos, sqrt, atan2, radians
import time

def find(parent, i):
    if parent[i] == i:
        return i
    return find(parent, parent[i])

def union(parent, rank, x, y):
    xroot = find(parent, x)
    yroot = find(parent, y)
    if rank[xroot] < rank[yroot]:
        parent[xroot] = yroot
    elif rank[xroot] > rank[yroot]:
        parent[yroot] = xroot
    else:
        parent[yroot] = xroot
        rank[xroot] += 1

def make_groups(indices):
    parent = {}
    rank = {}

    for i, j in indices:
        if i not in parent:
            parent[i] = i
            rank[i] = 0
        if j not in parent:
            parent[j] = j
            rank[j] = 0
        union(parent, rank, i, j)

    groups = {}
    for i in parent:
        root = find(parent, i)
        if root not in groups:
            groups[root] = []
        groups[root].append(i)

    return list(groups.values())

def calculate_distance(lat1, lon1, lat2, lon2):
    # Calculate the distance between two coordinates in kilometers
    R = 6371.0  # approximate radius of Earth in km
    dlon = radians(lon2) - radians(lon1)
    dlat = radians(lat2) - radians(lat1)
    a = sin(dlat / 2)**2 + cos(radians(lat1)) * cos(radians(lat2)) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    distance = R * c
    return distance

def remove_outliers(data, east_col='E.vel', north_col='N.vel', up_col='U.vel', threshold=0.2):
    # Calculate the magnitude of the velocity vectors
    magnitudes = np.sqrt(data[east_col] ** 2 + data[north_col] ** 2)

    # Calculate the azimuthal direction (in radians) of the velocity vectors
    azimuths = np.arctan2(data[north_col], data[east_col])

    # Calculate the median magnitude and median azimuth
    median_magnitude = np.median(magnitudes)
    median_azimuth = np.median(azimuths)

    # Calculate the magnitude and azimuthal differences from the median
    magnitude_diffs = np.abs(magnitudes - median_magnitude)
    azimuth_diffs = np.abs(np.arctan2(np.sin(azimuths - median_azimuth), np.cos(azimuths - median_azimuth)))

    # Compute the Interquartile Range (IQR) for both magnitude and azimuthal differences
    iqr_magnitude = np.percentile(magnitude_diffs, 75) - np.percentile(magnitude_diffs, 25)
    iqr_azimuth = np.percentile(azimuth_diffs, 75) - np.percentile(azimuth_diffs, 25)

    # Define the thresholds for outlier detection
    magnitude_threshold = threshold * iqr_magnitude
    azimuth_threshold = threshold * iqr_azimuth

    # Find the indices of stations with magnitude or azimuthal differences exceeding the thresholds
    outlier_indices = data.index[
        (magnitude_diffs > magnitude_threshold) |
        (azimuth_diffs > azimuth_threshold)
    ]

    # Check if all data points are outliers
    if len(outlier_indices) == len(data):
        # Return the median as there are no valid data points left
        median_velocities = data[[east_col, north_col, up_col]].median()
        return median_velocities, pd.DataFrame()

    # Remove the outliers from the dataset to get the data without outliers
    data_without_outliers = data.drop(outlier_indices)
    outliers = data.loc[outlier_indices]

    return data_without_outliers, outliers

def calculate_uncertainties(group_df):
    # Calculate the square root of the mean of squared uncertainties for each component
    squared_uncertainties = group_df[['E.sig', 'N.sig', 'U.sig']] ** 2
    mean_squared_uncertainties = squared_uncertainties.mean()
    uncertainties = np.sqrt(mean_squared_uncertainties)

    return uncertainties

def combine_velocities(input_folder, combined_folder):
    # Create the output folders if they don't exist
    os.makedirs(combined_folder, exist_ok=True)

    # Read all .vel files and merge them into a single velocity field
    file_paths = [f for f in os.listdir(input_folder) if f.endswith('.vel')]
    dfs = []
    for file_path in file_paths:
        basename = os.path.splitext(os.path.basename(file_path))[0]
        df = pd.read_csv(os.path.join(input_folder, file_path), delim_whitespace=True, header=None, skiprows=4)
        df.columns = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat']
        df['Ref'] = basename
        dfs.append(df)
    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df[(combined_df['Lon'] >= -15) & (combined_df['Lon'] <= 70) & (combined_df['Lat'] >= 5) & (combined_df['Lat'] <= 60)]

    # Create a separation matrix based on station distances
    stations = combined_df[['Lon', 'Lat']].values
    separation_matrix = pd.DataFrame(index=combined_df.index, columns=combined_df.index)
    for i, j in combinations(range(len(stations)), 2):
        distance = calculate_distance(stations[i][1], stations[i][0], stations[j][1], stations[j][0])
        separation_matrix.at[i, j] = separation_matrix.at[j, i] = distance

    # Identify stations within 0.1 km of each other and group them
    close_stations = separation_matrix[separation_matrix < 1.11].stack().index.tolist()
    close_stations_groups = make_groups(close_stations)
    
    #print(combined_df)
    
    for group in close_stations_groups:
        # Check if there is more than one station in the group
        if len(group) > 1:
            # Extract the relevant data for this group of stations
            group_df = combined_df.loc[group]
            group_df.columns = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat','Ref']
            #print("These are group_df contents:")
            #print(group_df)
            # Step 1: Remove outliers based on magnitude and azimuthal direction differences
            # For simplicity, we only consider the 'E.vel' and 'N.vel' components
            group_df[['E.vel', 'N.vel', 'U.vel']], outliers = remove_outliers(group_df[['E.vel', 'N.vel', 'U.vel']])

            # Step 2: Compute the median of the velocities
            median_velocities = group_df[['E.vel', 'N.vel', 'U.vel']].median().round(2)

            #print("These are group_df median velocities:")
            #print(median_velocities)
            
            # Step 3: Handle uncertainties
            # Calculate the square root of the sum of squared uncertainties for each component
            uncertainties = group_df[['E.sig', 'N.sig', 'U.sig']].median()
            uncertainties = uncertainties.round(2).astype('float')
            
            #print("These are group_df uncertainties:")
            #print(uncertainties)
            
            # Pick the first station in the group
            chosen_station_idx = 0
            chosen_station = group_df.iloc[chosen_station_idx]

            # Update the group_df with the combined values
            group_df[['E.vel', 'N.vel', 'U.vel']] = median_velocities
            group_df[['E.sig', 'N.sig', 'U.sig']] = uncertainties
            group_df['Lon'] = chosen_station['Lon'].round(5)
            group_df['Lat'] = chosen_station['Lat'].round(5)

            # Keep only the chosen station in the 'Stat' column
            group_df['Stat'] = chosen_station['Stat']
            
            # Assign 'E.adj', 'N.adj', 'U.adj', and 'Corr' values from the chosen station to the group DataFrame
            group_df['E.adj'] = chosen_station['E.adj'].round(2)
            group_df['N.adj'] = chosen_station['N.adj'].round(2)
            group_df['U.adj'] = chosen_station['U.adj'].round(2)
            group_df['Corr'] = chosen_station['Corr'].round(3)


            # Merge the processed group_df back into the combined_df
            combined_df.loc[group] = group_df
            
            # Additional debugging: Check if any NaN values exist in the merged DataFrame
            if combined_df.isnull().values.any():
                print("Warning: NaN values found in the merged DataFrame.")
                print(combined_df[combined_df.isnull().any(axis=1)])

        else:
            # If there is only one station in the group, just use it as is
            chosen_station = combined_df.loc[group[0]]

            # Keep only the chosen station in the 'Stat' column
            combined_df.loc[group, 'Stat'] = chosen_station['Stat']

    # Drop duplicates (keeping the first occurrence) from the combined_df based on 'Lon' and 'Lat'
    combined_df.drop_duplicates(subset=['Lon', 'Lat'], keep='first', inplace=True)
    
    # Drop the 'Ref' column from the combined dataframe
    combined_df.drop(columns=['Ref'], inplace=True)

    output_filename = "combined_vel_" + os.path.basename(input_folder)[-4:] + ".csv"
    combined_df.to_csv(os.path.join(combined_folder, output_filename), sep=' ', index=False)

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python combine_vel.py ./path2/input_folder ./path2/output_folder")
        sys.exit(1)

    input_folder = sys.argv[1]
    combined_folder = sys.argv[2] 

    # Time the execution of the combine_velocities function
    start_time = time.time()
    combine_velocities(input_folder, combined_folder)
    end_time = time.time()

    # Calculate and print the elapsed time
    elapsed_time = end_time - start_time
    print("Time taken to combine each data set: {:.2f} seconds".format(elapsed_time))
