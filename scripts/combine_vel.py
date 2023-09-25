import os
import sys
import pandas as pd
from itertools import product
import numpy as np
from math import sin, cos, sqrt, atan2, radians
import time
import warnings

# Ignore future warnings (I will fix these in a further update)
warnings.simplefilter(action='ignore', category=FutureWarning) 

# Implement a version of the Union-Find (also known as Disjoint Set) data structure. 
# The purpose of these functions is to track and merge groups of nearby GNSS stations

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

#def make_groups(indices):
#    parent = {}
#    rank = {}
#
#    for i, j in indices:
#        if i not in parent:
#            parent[i] = i
#            rank[i] = 0
#        if j not in parent:
#            parent[j] = j
#            rank[j] = 0
#        union(parent, rank, i, j)
#
#    groups = {}
#    for i in parent:
#        root = find(parent, i)
#        if root not in groups:
#            groups[root] = []
#        groups[root].append(i)
#
#    return list(groups.values())

# Modified version of make_groups function to include unconnected stations
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

    # Extract all unique indices from indices list
    all_indices = set(i for i, _ in indices) | set(j for _, j in indices)
    
    # Add stations that are not close to any other station.
    for idx in all_indices:
        if idx not in parent:
            groups[idx] = [idx]

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

def remove_outliers(data, east_col='E.vel', north_col='N.vel', up_col='U.vel'):
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
        # Return the median as there are no valid data points left
        median_velocities = data[[east_col, north_col, up_col]].median()
        return median_velocities, pd.DataFrame()

    # Remove the outliers from the dataset to get the data without outliers
    data_without_outliers = data.drop(outlier_indices)
    outliers = data.loc[outlier_indices]

    return data_without_outliers, outliers

# Testing faster approach to compute distances between stations
# The approach consists of creating a dictionary of stations that are within a certain distance of each other
# Insstead of computing the distance between all pairs of stations and saving them in a separation matrix

def create_distance_dict(stations, threshold=1.11):
    distance_dict = {}
    for i, j in product(range(len(stations)), repeat=2):
        distance = calculate_distance(stations[i][1], stations[i][0], stations[j][1], stations[j][0])
        if distance < threshold:
            if i not in distance_dict:
                distance_dict[i] = set()
            distance_dict[i].add(j)
    return distance_dict

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
    
    # Use the distance dictionary instead of a separation matrix
    distance_dict = create_distance_dict(stations)
    close_stations = [(i, j) for i, neighbours in distance_dict.items() for j in neighbours]
    close_stations_groups = make_groups(close_stations)
    
    #separation_matrix = pd.DataFrame(index=combined_df.index, columns=combined_df.index)

    #for i, j in product(range(len(stations)), repeat=2):
    #    distance = calculate_distance(stations[i][1], stations[i][0], stations[j][1], stations[j][0])
    #    separation_matrix.at[i, j] = distance

    #for i, j in combinations(range(len(stations)), 2):
    #    distance = calculate_distance(stations[i][1], stations[i][0], stations[j][1], stations[j][0])
    #    separation_matrix.at[i, j] = separation_matrix.at[j, i] = distance

    # Identify stations within 0.1 km of each other and group them
    #close_stations = separation_matrix[separation_matrix < 1.11].stack().index.tolist()
    #close_stations_groups = make_groups(close_stations)
    
    # Print the groups of close stations for debugging purposes
    #for index, group in enumerate(close_stations_groups, 1):
    #    print(f"Group {index} (Length: {len(group)}):\n {group}\n")

    #print(combined_df)
    # Create a folder called statistics inside the combined folder path to store the statistics of the combined velocity fields
    statistics_folder = os.path.join(combined_folder, "statistics")
    os.makedirs(statistics_folder, exist_ok=True)

    # Create a DataFrame to store the combined velocity fields
    aggregated_df = pd.DataFrame()

    # Create a DataFrame to store statistics of the combined velocity fields
    statistics_df = pd.DataFrame(columns=['Lon', 'Lat', 'Stat', 'Num'])
    
    for group in close_stations_groups:
        # Check if there is more than one station in the group
        if len(group) > 1:
            # Extract the relevant data for this group of stations
            group_df = combined_df.loc[group]
            group_df.columns = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat','Ref']
            #print("These are group_df contents:")
            #print(group_df)

            # Save the group_df to the aggregated_df DataFrame to be exported later as a CSV file for debugging purposes
            aggregated_df = pd.concat([aggregated_df, group_df], ignore_index=True)

            # Save group_df to a temorary variable
            # temp_df = group_df.copy()
            
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

            # Print df contents, median velocities, uncertainties and outliers in one print statement

            #print(f"Group:\n{temp_df}\nMedian velocities: {median_velocities}\nUncertainties: {uncertainties}\nOutliers: {outliers}\n")
            
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

            # Assign 'Ref' values from the chosen station to the group DataFrame
            # group_df['Ref'] = 'Clean'  # Set the 'Ref' column as 'Clean' for the chosen station

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
            # Print the group length and the group itself for debugging purposes
            # print(f"Group length: {len(group)}\nGroup:\n {group}")

            # If there is only one station in the group, just use it as is
            chosen_station = combined_df.loc[group[0]]

            # Keep only the chosen station in the 'Stat' column
            combined_df.loc[group, 'Stat'] = chosen_station['Stat']

            statistics_to_add = pd.DataFrame({
                'Lon': [chosen_station['Lon'].round(5)],
                'Lat': [chosen_station['Lat'].round(5)],
                'Stat': [chosen_station['Stat']],
                'Num': [1]
            })
            statistics_df = pd.concat([statistics_df, statistics_to_add], ignore_index=True)

            # Set the 'Ref' column as 'Clean'
            # combined_df.loc[group, 'Ref'] = 'Clean'

    # Drop duplicates (keeping the first occurrence) from the combined_df based on 'Lon' and 'Lat'
    combined_df.drop_duplicates(subset=['Lon', 'Lat'], keep='first', inplace=True)
    
    # Filter the dataframe to keep only the rows with 'Ref' as 'Clean'
    #combined_df = combined_df[combined_df['Ref'] == 'Clean']

    # Drop duplicates (keeping the first occurrence) from the combined_df based on 'Lon' and 'Lat'
    #combined_df = combined_df.drop_duplicates(subset=['Lat', 'Lon'])

    # Drop the 'Ref' column from the combined dataframe
    combined_df.drop(columns=['Ref'], inplace=True)

    # Save the combined_df to a CSV file
    output_filename = "combined_vel_" + os.path.basename(input_folder)[-4:] + ".csv"
    combined_df.to_csv(os.path.join(combined_folder, output_filename), sep=' ', index=False)

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
