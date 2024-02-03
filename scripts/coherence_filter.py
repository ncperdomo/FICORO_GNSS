import argparse
import json
import os
import sys
import glob
import pandas as pd
import numpy as np
from scipy.stats import lognorm
import concurrent.futures
import time

print() # Print a newline for better readability
print(f"################### Removing outliers using the Z-Score method ###################")

def haversine_distance(lon1, lat1, lon2, lat2):
    # Convert latitude and longitude from degrees to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    c = 2 * np.arcsin(np.sqrt(a))
    radius = 6371  # Earth's radius in kilometers
    return radius * c

def get_region_stringency(lon, lat, regions):
    default_sigma = 2  # Default sigma level
    for region in regions:
        if (region['min_lon'] <= lon <= region['max_lon'] and
                region['min_lat'] <= lat <= region['max_lat']):
            return region['sigma']
    return default_sigma

def filter_gps_velocities(file_name, radius=20, geo_strict=False, regions=[]):
    # Read the CSV file as a data frame, skipping the header row
    df = pd.read_csv(file_name, sep=' ', skiprows=1, header=None)
    df.columns = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat']
    
    # Empty set to store filtered stations
    filtered_stations = set()

    # Iterate over each GPS site
    for i, row in df.iterrows():
        site_lon, site_lat = row['Lon'], row['Lat']
        
        # Apply variable stringency if enabled, otherwise use default sigma level (2)
        sigma_level = get_region_stringency(site_lon, site_lat, regions) if geo_strict else 2

        # Calculate the Haversine distance between the GPS site and all stations
        distances = haversine_distance(site_lon, site_lat, df['Lon'].values, df['Lat'].values)
        # Filter stations that fall within the specified radius (strict adherence)
        nearby_stations = df[distances <= radius]

        # Proceed if there are more than 5 nearby stations
        if len(nearby_stations) >= 5:
            # Calculate mean and standard deviation of nearby stations' E.vel and N.vel
            e_vel_mean = nearby_stations['E.vel'].mean()
            e_vel_std = nearby_stations['E.vel'].std()
            n_vel_mean = nearby_stations['N.vel'].mean()
            n_vel_std = nearby_stations['N.vel'].std()

            e_vel_threshold = sigma_level * e_vel_std
            n_vel_threshold = sigma_level * n_vel_std
            # Filter stations with velocities outside the threshold
            filtered_stations.update(nearby_stations[
                (nearby_stations['E.vel'] < e_vel_mean - e_vel_threshold) |
                (nearby_stations['E.vel'] > e_vel_mean + e_vel_threshold) |
                (nearby_stations['N.vel'] < n_vel_mean - n_vel_threshold) |
                (nearby_stations['N.vel'] > n_vel_mean + n_vel_threshold)
            ].index)

    # Output results
    output_folder = './results/sites_excluded_coherence'
    os.makedirs(output_folder, exist_ok=True)
    output_clean_coherence = './results/output_coherence_analysis'
    os.makedirs(output_clean_coherence, exist_ok=True)

    # Save excluded stations
    filtered_df = df.loc[list(filtered_stations)].drop_duplicates()
    removed_lines_file = os.path.join(output_folder, f'{os.path.splitext(os.path.basename(file_name))[0]}.csv')
    filtered_df.to_csv(removed_lines_file, sep=' ', index=False)

    # Save included stations
    included_lines_df = df.drop(list(filtered_stations)).drop_duplicates()
    included_lines_file = os.path.join(output_clean_coherence, f'{os.path.splitext(os.path.basename(file_name))[0]}.csv')
    included_lines_df.to_csv(included_lines_file, sep=' ', index=False)

    # Printing results
    num_removed = len(filtered_df)
    num_total = len(df)
    percentage_removed = (num_removed / num_total) * 100
    text = f"\n----------------------------------------------------------------------------------\nNumber of stations removed for {os.path.basename(file_name)}: {num_removed} / {num_total} ({percentage_removed:.2f}%)\nSites excluded: {removed_lines_file}\nFiltered velocities: {included_lines_file}"
    print(text)

def parallel_filter_gps_velocities(folder_path, radius=20, geo_strict=False, regions=[]):
    # Find all CSV files in the folder
    file_names = glob.glob(os.path.join(folder_path, '*.csv'))
    # Create a ThreadPoolExecutor with the maximum number of worker threads
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Submit the filtering tasks for each file to the executor
        results = [executor.submit(filter_gps_velocities, file_name, radius, geo_strict, regions) for file_name in file_names]
        concurrent.futures.wait(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process GNSS data with optional geographic stringency.')
    parser.add_argument('folder_path', help='Path to the input folder containing CSV files')
    parser.add_argument('--geo_strict', action='store_true', help='Enable geographic stringency levels based on regions defined in a JSON file')
    parser.add_argument('--regions_json', type=str, help='Path to the JSON file with region definitions', default='')

    args = parser.parse_args()

    regions = []
    if args.geo_strict and args.regions_json:
        try:
            with open(args.regions_json, 'r') as file:
                regions = json.load(file)
        except FileNotFoundError:
            print(f"Error: JSON file {args.regions_json} not found.")
            sys.exit(1)
    
    # Time the execution of the parallel_filter_gps_velocities function
    start_time = time.time()
    parallel_filter_gps_velocities(args.folder_path, geo_strict=args.geo_strict, regions=regions)
    end_time = time.time()
    print(f"----------------------------------------------------------------------------------")
    print(f"Time taken: {end_time - start_time:.2f} seconds")
