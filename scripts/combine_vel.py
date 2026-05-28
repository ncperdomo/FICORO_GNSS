""" This code combines GNSS velocity fields into a single velocity field.
The input folder should contain a set of .vel files, previously cleaned
using the lognorm_filter and coherence_filter scripts. The output folder
will contain a single .csv file with the combined velocity field for each
reference frame. The code also outputs CSV files with the velocities
considered in the combination for each station and a file with the number
of solutions per GNSS station."""

import os
import sys
import pandas as pd
import numpy as np
from math import sin, cos, sqrt, atan2, radians
from sklearn.neighbors import BallTree
import time
import warnings

# Ignore future warnings (I will fix these in a future release)
warnings.simplefilter(action='ignore', category=FutureWarning)

""" Implement a version of the Union-Find (also known as Disjoint Set) data
structure. The purpose of these functions is to track and merge groups of
nearby GNSS stations"""

def find(parent, i):
    if parent[i] != i:
        parent[i] = find(parent, parent[i])  # path compression
    return parent[i]

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
    """ Build a neighbor dictionary mapping each station index to the set of
    station indices within `threshold` km. Uses a BallTree with the haversine
    metric for O(N log N) queries, replacing the original O(N²) pair loop."""
    coords_rad = np.radians(stations[:, [1, 0]])  # BallTree expects [lat, lon] in radians
    threshold_rad = threshold / 6371.0
    tree = BallTree(coords_rad, metric='haversine')
    indices = tree.query_radius(coords_rad, r=threshold_rad)
    return {i: set(nbrs.tolist()) for i, nbrs in enumerate(indices)}

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
        if basename.endswith('igb14'):
            df = pd.read_csv(os.path.join(input_folder, file_path), sep=r"\s+", header=None, skiprows=0)
        else:
            df = pd.read_csv(os.path.join(input_folder, file_path), sep=r"\s+", header=None, skiprows=4)
        df.columns = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat']
        df['Ref'] = basename
        dfs.append(df)
    combined_df = pd.concat(dfs, ignore_index=True)

    # Get the coordinates of all stations in the combined velocity field
    stations = combined_df[['Lon', 'Lat']].values

    # Use the distance dictionary to find close station pairs
    distance_dict = create_distance_dict(stations)
    close_stations = [(i, j) for i, neighbours in distance_dict.items() for j in neighbours]

    # Group close stations together
    close_stations_groups = make_groups(close_stations)
    print("Number of groups of close stations: {}".format(len(close_stations_groups)))

    # Create statistics folder
    statistics_folder = os.path.join(combined_folder, "statistics")
    os.makedirs(statistics_folder, exist_ok=True)

    # Column order for output
    _cols = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr',
             'U.vel', 'U.adj', 'U.sig', 'Stat']

    # Pre-extract numerical data and string arrays to avoid per-group pandas overhead.
    # Column layout: Lon=0, Lat=1, E.vel=2, N.vel=3, E.adj=4, N.adj=5,
    #                E.sig=6, N.sig=7, Corr=8, U.vel=9, U.adj=10, U.sig=11
    _num_cols = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj',
                 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig']
    data_np = combined_df[_num_cols].to_numpy(dtype=np.float64)
    stat_arr = combined_df['Stat'].to_numpy()
    ref_arr  = combined_df['Ref'].to_numpy()

    # Accumulate results as lists; build DataFrames once after the loop
    output_rows = []
    statistics_rows = []
    aggregated_frames = []

    for group in close_stations_groups:
        g = np.asarray(group, dtype=np.intp)
        rows = data_np[g]          # (k, 12) numpy slice — no pandas overhead
        chosen_row = rows[0]
        chosen_ref  = ref_arr[g[0]]
        chosen_stat = stat_arr[g[0]]

        if len(g) > 1:
            # Save pre-outlier-removal data for eura debugging
            if chosen_ref.endswith('eura'):
                aggregated_frames.append(combined_df.iloc[g])

            # --- Outlier detection (IQR on magnitude/azimuth) ---
            # In the normal case (not all-outliers), the original code's
            # group_df.loc[result.index] = result.values is a no-op, so the
            # medians below are always computed over the full un-modified rows.
            # We therefore compute medians directly and only invoke the outlier
            # test to emit the same warning messages as the original code.
            e_vel = rows[:, 2]
            n_vel = rows[:, 3]
            u_vel = rows[:, 9]

            magnitudes    = np.sqrt(e_vel**2 + n_vel**2)
            azimuths      = np.arctan2(n_vel, e_vel)
            med_mag       = np.median(magnitudes)
            med_az        = np.median(azimuths)
            mag_diffs     = np.abs(magnitudes - med_mag)
            az_diffs      = np.abs(np.arctan2(np.sin(azimuths - med_az),
                                               np.cos(azimuths - med_az)))
            q1_m, q3_m    = np.percentile(mag_diffs, [25, 75])
            q1_a, q3_a    = np.percentile(az_diffs,  [25, 75])
            iqr_m, iqr_a  = q3_m - q1_m, q3_a - q1_a
            is_outlier    = (
                (mag_diffs < q1_m - 1.5*iqr_m) | (mag_diffs > q3_m + 1.5*iqr_m) |
                (az_diffs  < q1_a - 1.5*iqr_a) | (az_diffs  > q3_a + 1.5*iqr_a)
            )

            if is_outlier.all():
                # All-outliers edge case: warn if all U values are also zero
                if not np.any(u_vel != 0):
                    print("Warning: All vertical velocities are zero. Assigning 0.00 as the median.")

            # Step 2: Median velocities (over all rows — see note above)
            # Use np.round (not Python built-in round) to match pandas .round() behavior
            e_med = float(np.round(np.median(e_vel), 2))
            n_med = float(np.round(np.median(n_vel), 2))
            nonzero_u = u_vel[u_vel != 0]
            if len(nonzero_u) == 0:
                print("Warning: All vertical velocities are zero. Assigning NaN as the median.")
                u_median = np.nan
            else:
                u_median = float(np.round(np.median(nonzero_u), 2))

            # Step 3: Median uncertainties
            e_sig = float(np.round(np.median(rows[:, 6]), 2))
            n_sig = float(np.round(np.median(rows[:, 7]), 2))
            u_sig = float(np.round(np.median(rows[:, 11]), 2))

            output_rows.append({
                'Lon':   float(np.round(chosen_row[0], 5)),
                'Lat':   float(np.round(chosen_row[1], 5)),
                'E.vel': e_med,
                'N.vel': n_med,
                'E.adj': float(np.round(chosen_row[4], 2)),
                'N.adj': float(np.round(chosen_row[5], 2)),
                'E.sig': e_sig,
                'N.sig': n_sig,
                'Corr':  float(np.round(chosen_row[8], 3)),
                'U.vel': u_median,
                'U.adj': float(np.round(chosen_row[10], 2)),
                'U.sig': u_sig,
                'Stat':  chosen_stat,
            })

            # Record statistics for eura reference frame
            if chosen_ref.endswith('eura'):
                statistics_rows.append({
                    'Lon':  float(np.round(chosen_row[0], 5)),
                    'Lat':  float(np.round(chosen_row[1], 5)),
                    'Stat': chosen_stat,
                    'Num':  len(g),
                })

        else:
            # Single station — use as is
            output_rows.append({
                'Lon':   float(chosen_row[0]),
                'Lat':   float(chosen_row[1]),
                'E.vel': float(chosen_row[2]),
                'N.vel': float(chosen_row[3]),
                'E.adj': float(chosen_row[4]),
                'N.adj': float(chosen_row[5]),
                'E.sig': float(chosen_row[6]),
                'N.sig': float(chosen_row[7]),
                'Corr':  float(chosen_row[8]),
                'U.vel': float(chosen_row[9]),
                'U.adj': float(chosen_row[10]),
                'U.sig': float(chosen_row[11]),
                'Stat':  chosen_stat,
            })

            if chosen_ref.endswith('eura'):
                statistics_rows.append({
                    'Lon':  float(np.round(chosen_row[0], 5)),
                    'Lat':  float(np.round(chosen_row[1], 5)),
                    'Stat': chosen_stat,
                    'Num':  1,
                })

    # Build output DataFrame from accumulated rows (avoids O(N²) concat-in-loop)
    output_df = pd.DataFrame(output_rows, columns=_cols)

    # Debugging: check for unexpected NaN values (moved outside the group loop)
    if output_df.isnull().values.any():
        print("Warning: NaN values found in the merged DataFrame.")
        print(output_df[output_df.isnull().any(axis=1)])

    # Determine output filename
    if basename.endswith('igb14'):
        output_filename = "combined_vel_igb14.csv"
    else:
        output_filename = "combined_vel_" + os.path.basename(input_folder)[-4:] + ".csv"

    output_df.to_csv(os.path.join(combined_folder, output_filename), sep=' ', index=False)

    # Save debugging and statistics files if eura reference frame was processed
    if aggregated_frames:
        aggregated_df = pd.concat(aggregated_frames, ignore_index=True)
        aggregated_df.to_csv(os.path.join(statistics_folder, "grouped_stations.csv"), sep=',', index=False)

    if statistics_rows:
        statistics_df = pd.DataFrame(statistics_rows, columns=['Lon', 'Lat', 'Stat', 'Num'])
        statistics_df.to_csv(os.path.join(statistics_folder, "site_statistics.csv"), sep=',', index=False)

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
