import os
import sys 
import glob
import pandas as pd
import numpy as np
from scipy.stats import lognorm
import concurrent.futures
import time
import random

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
    distance = radius * c

    return distance

def filter_gps_velocities(file_name, radius=20):
    # Read the CSV file as a data frame, skipping the header row
    df = pd.read_csv(file_name, sep=' ', skiprows=1, header=None)
    df.columns = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat']

    # Empty set to store filtered stations
    filtered_stations = set()

    # Iterate over each GPS site
    for i, row in df.iterrows():
        site_lon, site_lat = row['Lon'], row['Lat']

        # Calculate the Haversine distance between the GPS site and all stations
        distances = haversine_distance(site_lon, site_lat, df['Lon'], df['Lat'])

        # Filter stations that fall within the specified radius (strict adherence)
        nearby_stations = df[distances <= radius]

        # Proceed if there are more than 5 nearby stations
        if len(nearby_stations) >= 5:
            # Calculate mean and standard deviation of nearby stations' E.vel and N.vel
            e_vel_mean = nearby_stations['E.vel'].mean()
            e_vel_std = nearby_stations['E.vel'].std()

            n_vel_mean = nearby_stations['N.vel'].mean()
            n_vel_std = nearby_stations['N.vel'].std()

            # Check if the CSV file name contains the word "hamiel"
            if "hamiel" in file_name:
                # Apply 3 standard deviations for "hamiel" files
                e_vel_threshold = 3 * e_vel_std
                n_vel_threshold = 3 * n_vel_std
            else:
                # Apply 2 standard deviations for other files
                e_vel_threshold = 2 * e_vel_std
                n_vel_threshold = 2 * n_vel_std

            # Filter stations with velocities outside the threshold
            filtered_stations.update(nearby_stations[
                (nearby_stations['E.vel'] < e_vel_mean - e_vel_threshold) |
                (nearby_stations['E.vel'] > e_vel_mean + e_vel_threshold) |
                (nearby_stations['N.vel'] < n_vel_mean - n_vel_threshold) |
                (nearby_stations['N.vel'] > n_vel_mean + n_vel_threshold)
            ].index)

    # Create a data frame with the filtered stations
    filtered_df = df.loc[list(filtered_stations)].drop_duplicates()

# Print the number of removed stations for the current dataset
    file_name = os.path.splitext(os.path.basename(file_name))[0]
    num_removed = len(filtered_df)
    num_total = len(df)
    percentage_removed = (num_removed / num_total) * 100

    #print(f"----------------------------------------------------------------------------------")
    #print(f"Number of stations removed for {file_name}: {num_removed} / {num_total} ({percentage_removed:.2f}%)")

    # Create a directory to store the files listing filtered stations
    output_folder = './results/sites_excluded_coherence'
    os.makedirs(output_folder, exist_ok=True)

    # Create a directory to store the clean files listing stations that were not removed
    output_clean_coherence = './results/output_coherence_analysis'
    os.makedirs(output_clean_coherence, exist_ok=True)

    # Create a data frame including the lines in the original df, removing those in filtered_df
    included_lines_df = df.drop(list(filtered_stations)).drop_duplicates()

    # Save the filtered stations to a new CSV file
    removed_lines_file = os.path.join(output_folder, f'{file_name}.csv')
    filtered_df.to_csv(removed_lines_file, sep=' ', index=False)
    #print(f"Sites excluded: {removed_lines_file}")

    # Save the included lines to a new CSV file
    included_lines_file = os.path.join(output_clean_coherence, f'{file_name}.csv')
    included_lines_df.to_csv(included_lines_file, sep=' ', index=False)
    #print(f"Filtered velocities: {included_lines_file}")

    # Define the text to be printed in multiple lines
    text = f"\n----------------------------------------------------------------------------------\nNumber of stations removed for {file_name}: {num_removed} / {num_total} ({percentage_removed:.2f}%)\nSites excluded: {removed_lines_file}\nFiltered velocities: {included_lines_file}"
    print(text)
    


def parallel_filter_gps_velocities(folder_path, radius=20):
    # Find all CSV files in the folder
    file_names = glob.glob(os.path.join(folder_path, '*.csv'))

    # Create a ThreadPoolExecutor with the maximum number of worker threads
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # Submit the filtering tasks for each file to the executor
        results = [executor.submit(filter_gps_velocities, file_name, radius) for file_name in file_names]

        # Wait for all tasks to complete
        concurrent.futures.wait(results)

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python coherence_filter.py ./path2/input_folder")
        sys.exit(1)

    folder_path = sys.argv[1]

    # Time the execution of the combine_velocities function
    start_time = time.time()
    parallel_filter_gps_velocities(folder_path)
    end_time = time.time()

    # Calculate and print the elapsed time
    elapsed_time = end_time - start_time
    print("Time taken to run coherence_filter.py: {:.2f} seconds".format(elapsed_time))
