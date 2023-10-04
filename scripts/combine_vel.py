""" This code combines GNSS velocity fields into a single velocity field. 
The input folder should contain a set of .vel files, previously cleaned 
using the lognorm_filter and coherence_filter scripts. The output folder 
will contain a single .csv file with the combined velocity field for each 
reference frame. The code also outputs CSV files with the velocities 
considered in the combination for each station and a file with the number 
of solutions per GNSS station."""

""" Import necessary modules """
import os
import sys
import pandas as pd
from itertools import product
import numpy as np
from math import sin, cos, sqrt, atan2, radians
import time
import warnings

# Ignore future warnings (I will fix these in a future release)
warnings.simplefilter(action='ignore', category=FutureWarning) 

""" Implement a version of the Union-Find (also known as Disjoint Set) data 
structure. The purpose of these functions is to track and merge groups of
nearby GNSS stations""" 

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

# Modified version of make_groups function to include unconnected stations

def make_groups(indices):
    """ make_groups is a function that uses the above Union-Find implementation 
    to group GNSS stations based on their proximity. It takes a list of indices
    as input and returns a list of groups of indices. Each group contains the 
    indices of stations that are close to each other."""

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

    # Extract all unique indices from indices list
    all_indices = set(i for i, _ in indices) | set(j for _, j in indices)
    
    # Add stations that are not close to any other station.
    for idx in all_indices:
        if idx not in parent:
            groups[idx] = [idx]

    return list(groups.values())

def calculate_distance(lat1, lon1, lat2, lon2):
    """ The calculate_distance function computes the Haversine distance between two sets 
    of latitude and longitude values, returning the result in kilometers."""

    # Calculate the distance between two coordinates in kilometers
    R = 6371.0  # approximate radius of Earth in km
    dlon = radians(lon2) - radians(lon1)
    dlat = radians(lat2) - radians(lat1)
    a = sin(dlat / 2)**2 + cos(radians(lat1)) * cos(radians(lat2)) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    distance = R * c
    return distance

def remove_outliers(data, east_col='E.vel', north_col='N.vel', up_col='U.vel'):
    """ The remove_outliers function removes outliers from the dataset based on the 
    magnitude and azimuthal direction of the velocity vectors. It takes a DataFrame 
    as input and returns the data without outliers and the outliers. The function
    implements the Interquartile Range (IQR) method to detect outliers."""

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
    Q1_magnitude = np.percentile(magnitude_diffs, 25)
    Q3_magnitude = np.percentile(magnitude_diffs, 75)
    iqr_magnitude = Q3_magnitude - Q1_magnitude

    Q1_azimuth = np.percentile(azimuth_diffs, 25)
    Q3_azimuth = np.percentile(azimuth_diffs, 75)
    iqr_azimuth = Q3_azimuth - Q1_azimuth

    # Define the thresholds for outlier detection
    lower_magnitude_threshold = Q1_magnitude - 1.5 * iqr_magnitude
    upper_magnitude_threshold = Q3_magnitude + 1.5 * iqr_magnitude

    lower_azimuth_threshold = Q1_azimuth - 1.5 * iqr_azimuth
    upper_azimuth_threshold = Q3_azimuth + 1.5 * iqr_azimuth

    # Find the indices of stations with magnitude or azimuthal differences exceeding the thresholds
    outlier_indices = data.index[
        (magnitude_diffs < lower_magnitude_threshold) |
        (magnitude_diffs > upper_magnitude_threshold) |
        (azimuth_diffs < lower_azimuth_threshold) |
        (azimuth_diffs > upper_azimuth_threshold)
    ]

    # Check if all data points are outliers
    if len(outlier_indices) == len(data):
        # Compute horizontal and vertical median velocities separately. 
        # Return the median East and North velocity components as there are no valid data points left
        median_velocities = data[[east_col, north_col]].median()
        # For the vertical component, consider only non-zero values in the median calculation
        median_velocities[up_col] = data[data[up_col] != 0][up_col].median()
        # If all values are zero, return 0.00 as the median (later, the code will detect zero values and assign NaN)
        if data[data[up_col] != 0][up_col].empty:
            print("Warning: All vertical velocities are zero. Assigning 0.00 as the median.")
            median_velocities[up_col] = round(0.00,2)
        return median_velocities, pd.DataFrame()

    # Remove the outliers from the dataset to get the data without outliers
    data_without_outliers = data.drop(outlier_indices)
    outliers = data.loc[outlier_indices]

    # Return the data without outliers and the outliers
    return data_without_outliers, outliers

def create_distance_dict(stations, threshold=1.11):
    """ Instead of creating a separation matrix for station distances, a dictionary 
    approach is used to efficiently map stations within a certain distance of each other.
    This approach reduces the time complexity of the algorithm from O(n^2) to O(n)."""
    distance_dict = {}
    for i, j in product(range(len(stations)), repeat=2):
        distance = calculate_distance(stations[i][1], stations[i][0], stations[j][1], stations[j][0])
        if distance < threshold:
            if i not in distance_dict:
                distance_dict[i] = set()
            distance_dict[i].add(j)
    return distance_dict

def combine_velocities(input_folder, combined_folder):
    """ The combine_velocities function takes an input folder path containing previously 
    filtered .vel files and an output folder path, where the combined velocity field in 
    different reference frames will be saved. The combination is done by:
    - Reading multiple .vel files and merging their data.
    - Creating a distance dictionary that maps station pairs based on their proximity.
    - Using the distance dictionary, it groups close stations together.
    For each group of close stations, it:
        - Removes outliers from the group based on magnitude and azimuthal direction differences.
        - Computes the median of the velocities and uncertainties for each component.
        - Updates the velocity and other fields for the group based on the first station in the group.
        - Records statistics for the group (number of solutions per station)
    - After processing all groups, it saves the combined velocity field as a .csv file"""

    # Create the output folders if they don't exist
    os.makedirs(combined_folder, exist_ok=True)

    # Read all .vel files and merge them into a single velocity field
    file_paths = [f for f in os.listdir(input_folder) if f.endswith('.vel')]
    dfs = []
    for file_path in file_paths:
        basename = os.path.splitext(os.path.basename(file_path))[0]
        # If basename ends with igb14 set skiprows to 0, otherwise set skiprows to 4, because the igb14 files have no header
        if basename.endswith('igb14'):
            df = pd.read_csv(os.path.join(input_folder, file_path), delim_whitespace=True, header=None, skiprows=0)
        else:
            df = pd.read_csv(os.path.join(input_folder, file_path), delim_whitespace=True, header=None, skiprows=4)
        df.columns = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat']
        df['Ref'] = basename
        dfs.append(df)
    combined_df = pd.concat(dfs, ignore_index=True)
    # combined_df[(combined_df['Lon'] >= -15) & (combined_df['Lon'] <= 70) & (combined_df['Lat'] >= 5) & (combined_df['Lat'] <= 60)]

    # Get the coordinates of all stations in the combined velocity field as a numpy array of shape (n, 2) where n is the number of stations 
    stations = combined_df[['Lon', 'Lat']].values
    
    # Use the distance dictionary instead of a separation matrix to reduce the time complexity of the algorithm
    distance_dict = create_distance_dict(stations)
    close_stations = [(i, j) for i, neighbours in distance_dict.items() for j in neighbours] # List of tuples of close station pairs

    # Group close stations together based on the distance dictionary
    close_stations_groups = make_groups(close_stations) # List of lists of close stations

    # Check the length of the close_stations_groups list
    print("Number of groups of close stations: {}".format(len(close_stations_groups)))

    # Create a folder called statistics inside the combined folder path to store the statistics of the combined velocity fields
    statistics_folder = os.path.join(combined_folder, "statistics")
    os.makedirs(statistics_folder, exist_ok=True)

    # Create a DataFrame to store the combined velocity fields. Only if basename ends with eura
    if basename.endswith('eura'):
        aggregated_df = pd.DataFrame()

    # Create a DataFrame to store statistics of the combined velocity field. Only if basename ends with eura
    if basename.endswith('eura'):
        statistics_df = pd.DataFrame(columns=['Lon', 'Lat', 'Stat', 'Num'])
    
    for group in close_stations_groups:
        # Check if there is more than one station in the group
        if len(group) > 1:
            # Extract the relevant data for this group of stations
            group_df = combined_df.loc[group]
            group_df.columns = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat','Ref']

            # Save the group_df to the aggregated_df DataFrame to be exported later as a CSV file for debugging purposes. Only if basename ends with eura
            if group_df['Ref'].iloc[0].endswith('eura'):
                aggregated_df = pd.concat([aggregated_df, group_df], ignore_index=True)
            
            # Step 1: Remove outliers based on magnitude and azimuthal direction differences
            # For simplicity, we only consider the 'E.vel' and 'N.vel' components
            group_df[['E.vel', 'N.vel', 'U.vel']], outliers = remove_outliers(group_df[['E.vel', 'N.vel', 'U.vel']])

            # Step 2: Compute the median of horzontal and vertical velocities separately
            # For the vertical component, we only include non-zero values in the median calculation. 
            median_velocities = group_df[['E.vel', 'N.vel']].median().round(2)
            median_velocities['U.vel'] = group_df[group_df['U.vel'] != 0]['U.vel'].median().round(2)
            # If all vertical velocities are zero (i.e., the input velocity fields did not estimate verticals), return NaN as the median.
            if group_df[group_df['U.vel'] != 0]['U.vel'].empty:
                print("Warning: All vertical velocities are zero. Assigning NaN as the median.")
                median_velocities['U.vel'] = np.nan
            
            # Step 3: Compute median uncertainties for each velocity component
            uncertainties = group_df[['E.sig', 'N.sig', 'U.sig']].median()
            uncertainties = uncertainties.round(2).astype('float')
            
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

            # Save the number of stations in the group to the statistics_df DataFrame if basename ends with eura
            if chosen_station['Ref'].endswith('eura'):
                statistics_to_add = pd.DataFrame({
                    'Lon': [chosen_station['Lon'].round(5)],
                    'Lat': [chosen_station['Lat'].round(5)],
                    'Stat': [chosen_station['Stat']],
                    'Num': [len(group)]
                })
                statistics_df = pd.concat([statistics_df, statistics_to_add], ignore_index=True)

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

            # Save the number of stations in the group [1] to the statistics_df DataFrame if basename ends with eura
            if chosen_station['Ref'].endswith('eura'):
                statistics_to_add = pd.DataFrame({
                    'Lon': [chosen_station['Lon'].round(5)],
                    'Lat': [chosen_station['Lat'].round(5)],
                    'Stat': [chosen_station['Stat']],
                    'Num': [1]
                })
                statistics_df = pd.concat([statistics_df, statistics_to_add], ignore_index=True)

    # Drop duplicates (keeping the first occurrence) from the combined_df based on 'Lon' and 'Lat'
    combined_df.drop_duplicates(subset=['Lon', 'Lat'], keep='first', inplace=True)

    # Drop the 'Ref' column from the combined dataframe
    combined_df.drop(columns=['Ref'], inplace=True)

    # Save the combined_df to a CSV file
    # If basename ends with igb14, set the output filename to combined_vel_igb14.csv, 
    # otherwise set based on the last 4 characters of the input folder name
    if basename.endswith('igb14'):
        output_filename = "combined_vel_igb14.csv"
    else:
        output_filename = "combined_vel_" + os.path.basename(input_folder)[-4:] + ".csv"
    
    # Save the combined velocity field to a CSV file
    combined_df.to_csv(os.path.join(combined_folder, output_filename), sep=' ', index=False)

    if chosen_station['Ref'].endswith('eura'):
        # Save groupped stations to a CSV file for debugging purposes
        group_df_file_path = os.path.join(statistics_folder, "grouped_stations.csv")
        aggregated_df.to_csv(group_df_file_path, sep=',', index=False)

        # Save the statistics_df to a CSV file
        statistics_df_file_path = os.path.join(statistics_folder, "site_statistics.csv")
        statistics_df.to_csv(statistics_df_file_path, sep=',', index=False)

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

    # Calculate and print the elapsed time in minutes
    elapsed_time = (end_time - start_time) / 60
    print("Time taken to combine GNSS velocity fields: {:.2f} minutes".format(elapsed_time))
