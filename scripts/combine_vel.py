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
from sklearn.neighbors import BallTree
import time

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


def create_distance_dict(stations, threshold=1.11):
    """Build a neighbor dictionary mapping each station index to the set of
    station indices within `threshold` km. Uses a BallTree with the haversine
    metric for O(N log N) queries, replacing the original O(N^2) pair loop."""
    coords_rad = np.radians(stations[:, [1, 0]])  # BallTree expects [lat, lon] in radians
    threshold_rad = threshold / 6371.0
    tree = BallTree(coords_rad, metric='haversine')
    indices = tree.query_radius(coords_rad, r=threshold_rad)
    return {i: set(nbrs.tolist()) for i, nbrs in enumerate(indices)}


def combine_component(values, sigmas, combine_method="median", ignore_zeros=False):
    """Combine a single component using either the median or an inverse-variance weighted mean.

    Parameters
    ----------
    values : array-like
        Component values to combine.
    sigmas : array-like
        Uncertainty values corresponding to `values`.
    combine_method : str
        Either 'median' or 'weighted_mean'.
    ignore_zeros : bool
        If True, drop entries where the component value is exactly zero before combining.
        This is used for the vertical velocity component.

    Returns
    -------
    tuple[float, float]
        Combined value and combined uncertainty.
    """
    values = np.asarray(values, dtype=np.float64)
    sigmas = np.asarray(sigmas, dtype=np.float64)

    valid = np.isfinite(values) & np.isfinite(sigmas)
    values = values[valid]
    sigmas = sigmas[valid]

    if ignore_zeros:
        nonzero = values != 0
        values = values[nonzero]
        sigmas = sigmas[nonzero]

    if len(values) == 0:
        return np.nan, np.nan

    if combine_method == "median":
        return float(np.round(np.median(values), 2)), float(np.round(np.median(sigmas), 2))

    if combine_method == "weighted_mean":
        valid_sigma = np.isfinite(sigmas) & (sigmas > 0)
        values = values[valid_sigma]
        sigmas = sigmas[valid_sigma]

        if len(values) == 0:
            return np.nan, np.nan

        weights = 1.0 / np.square(sigmas)
        weight_sum = np.sum(weights)

        if weight_sum == 0:
            return np.nan, np.nan

        combined_value = np.sum(weights * values) / weight_sum
        combined_sigma = np.sqrt(1.0 / weight_sum)
        return float(np.round(combined_value, 2)), float(np.round(combined_sigma, 2))

    raise ValueError("combine_method must be either 'median' or 'weighted_mean'")


def combine_velocities(input_folder, combined_folder, combine_method="median"):
    """Combine previously filtered .vel files into a single velocity field.

    Parameters
    ----------
    input_folder : str
        Folder containing the input .vel files.
    combined_folder : str
        Folder where the combined velocity field and statistics will be saved.
    combine_method : str
        Either 'median' or 'weighted_mean'.
    """

    if combine_method not in {"median", "weighted_mean"}:
        raise ValueError("combine_method must be either 'median' or 'weighted_mean'")

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
    _num_cols = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj',
                 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig']
    data_np = combined_df[_num_cols].to_numpy(dtype=np.float64)
    stat_arr = combined_df['Stat'].to_numpy()
    ref_arr = combined_df['Ref'].to_numpy()

    # Accumulate results as lists; build DataFrames once after the loop
    output_rows = []
    statistics_rows = []
    aggregated_frames = []

    for group in close_stations_groups:
        g = np.asarray(group, dtype=np.intp)
        rows = data_np[g]          # (k, 12) numpy slice — no pandas overhead
        chosen_row = rows[0]
        chosen_ref = ref_arr[g[0]]
        chosen_stat = stat_arr[g[0]]

        if len(g) > 1:
            # Save pre-combination data for eura debugging
            if chosen_ref.endswith('eura'):
                aggregated_frames.append(combined_df.iloc[g])

            # Combine the group using the selected method
            e_vel, e_sig = combine_component(rows[:, 2], rows[:, 6], combine_method=combine_method)
            n_vel, n_sig = combine_component(rows[:, 3], rows[:, 7], combine_method=combine_method)
            u_vel, u_sig = combine_component(rows[:, 9], rows[:, 11], combine_method=combine_method, ignore_zeros=True)

            # Spread diagnostic: max Euclidean distance (in velocity units) of
            # any solution from the group's combined horizontal velocity. This
            # flags groups whose co-located solutions disagree badly (likely
            # offsets, mislabels, or two monuments inside the merge radius) for
            # later inspection, without affecting the combined velocity.
            spread = float(np.round(
                np.max(np.sqrt((rows[:, 2] - e_vel)**2 +
                               (rows[:, 3] - n_vel)**2)), 2))

            output_rows.append({
                'Lon': float(np.round(chosen_row[0], 5)),
                'Lat': float(np.round(chosen_row[1], 5)),
                'E.vel': e_vel,
                'N.vel': n_vel,
                'E.adj': float(np.round(chosen_row[4], 2)),
                'N.adj': float(np.round(chosen_row[5], 2)),
                'E.sig': e_sig,
                'N.sig': n_sig,
                'Corr': float(np.round(chosen_row[8], 3)),
                'U.vel': u_vel,
                'U.adj': float(np.round(chosen_row[10], 2)),
                'U.sig': u_sig,
                'Stat': chosen_stat,
            })

            # Record statistics for eura reference frame
            if chosen_ref.endswith('eura'):
                statistics_rows.append({
                    'Lon': float(np.round(chosen_row[0], 5)),
                    'Lat': float(np.round(chosen_row[1], 5)),
                    'Stat': chosen_stat,
                    'Num': len(g),
                    'Spread': spread,
                })

        else:
            # Single station — use as is
            output_rows.append({
                'Lon': float(chosen_row[0]),
                'Lat': float(chosen_row[1]),
                'E.vel': float(chosen_row[2]),
                'N.vel': float(chosen_row[3]),
                'E.adj': float(chosen_row[4]),
                'N.adj': float(chosen_row[5]),
                'E.sig': float(chosen_row[6]),
                'N.sig': float(chosen_row[7]),
                'Corr': float(chosen_row[8]),
                'U.vel': float(chosen_row[9]),
                'U.adj': float(chosen_row[10]),
                'U.sig': float(chosen_row[11]),
                'Stat': chosen_stat,
            })

            if chosen_ref.endswith('eura'):
                statistics_rows.append({
                    'Lon': float(np.round(chosen_row[0], 5)),
                    'Lat': float(np.round(chosen_row[1], 5)),
                    'Stat': chosen_stat,
                    'Num': 1,
                    'Spread': 0.0,
                })

    # Build output DataFrame from accumulated rows (avoids O(N^2) concat-in-loop)
    output_df = pd.DataFrame(output_rows, columns=_cols)

    # Debugging: check for unexpected NaN values (moved outside the group loop)
    if output_df.isnull().values.any():
        print("Warning: NaN values found in the merged DataFrame.")
        print(output_df[output_df.isnull().any(axis=1)])

    # Determine output filename from the run-level folder name
    folder_name = os.path.basename(os.path.normpath(input_folder))
    if folder_name.endswith('igb14'):
        output_filename = "combined_vel_igb14.csv"
    else:
        output_filename = "combined_vel_" + folder_name[-4:] + ".csv"

    output_df.to_csv(os.path.join(combined_folder, output_filename), sep=' ', index=False)

    # Save debugging and statistics files if eura reference frame was processed
    if aggregated_frames:
        aggregated_df = pd.concat(aggregated_frames, ignore_index=True)
        aggregated_df.to_csv(os.path.join(statistics_folder, "grouped_stations.csv"), sep=',', index=False)

    if statistics_rows:
        statistics_df = pd.DataFrame(statistics_rows, columns=['Lon', 'Lat', 'Stat', 'Num', 'Spread'])
        statistics_df.to_csv(os.path.join(statistics_folder, "site_statistics.csv"), sep=',', index=False)


if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) not in (3, 4):
        print("Usage: python combine_vel.py ./path2/input_folder ./path2/output_folder [median|weighted_mean]")
        sys.exit(1)

    input_folder = sys.argv[1]
    combined_folder = sys.argv[2]
    combine_method = sys.argv[3] if len(sys.argv) == 4 else "median"

    if combine_method not in {"median", "weighted_mean"}:
        print("Error: combine_method must be either 'median' or 'weighted_mean'")
        sys.exit(1)

    # Time the execution of the combine_velocities function
    start_time = time.time()
    combine_velocities(input_folder, combined_folder, combine_method=combine_method)
    end_time = time.time()

    # Calculate and print the elapsed time in minutes
    elapsed_time = (end_time - start_time)
    print("Time taken to combine GNSS velocity fields: {:.2f} seconds".format(elapsed_time))