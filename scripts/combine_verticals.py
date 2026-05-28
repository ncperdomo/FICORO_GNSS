""" Combines vertical GNSS velocity fields from multiple .vel files into a single
combined field. Input files must have 5 space-separated columns (no header):
    Lon  Lat  U.vel  U.sig  Stat
Output: combined_vertical_velocity_field.vel (same 5-column format, no header)
        debug_combined_vertical_velocity_field.log
"""

import os
import sys
import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree
import time


# ---------------------------------------------------------------------------
# Union-Find helpers (path-compressed)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Spatial neighbor search
# ---------------------------------------------------------------------------

def create_distance_dict(stations, threshold=1.2):
    """Return a dict mapping each station index to the set of neighbor indices
    within `threshold` km, excluding self-pairs.

    Uses a BallTree with the haversine metric — O(N log N) — replacing the
    original O(N²) itertools.product loop.

    Longitudes are normalized to [-180, 180] before the search so that input
    files using the 0–360° convention are correctly matched against files using
    the -180–180° convention (e.g. serpelloni_2022 / castro_2021)."""
    # Normalize longitudes to [-180, 180] for consistent haversine distances
    lons_norm = ((stations[:, 0] + 180.0) % 360.0) - 180.0
    coords_rad = np.radians(np.column_stack([stations[:, 1], lons_norm]))  # [lat, lon]
    threshold_rad = threshold / 6371.0
    tree = BallTree(coords_rad, metric='haversine')
    raw = tree.query_radius(coords_rad, r=threshold_rad)
    return {
        i: set(nbrs.tolist()) - {i}
        for i, nbrs in enumerate(raw)
        if len(nbrs) > 1  # len > 1 means there is at least one neighbor besides self
    }


# ---------------------------------------------------------------------------
# Union-Find grouping
# ---------------------------------------------------------------------------

def make_groups(distance_dict, n):
    """Group station indices using the neighbor dict and Union-Find.
    Isolated stations (no neighbors within threshold) are returned as
    singleton groups."""
    parent = {i: i for i in range(n)}
    rank   = {i: 0 for i in range(n)}

    for i, neighbors in distance_dict.items():
        for j in neighbors:
            union(parent, rank, i, j)

    groups = {}
    for i in parent:
        root = find(parent, i)
        if root not in groups:
            groups[root] = []
        groups[root].append(i)

    # Add stations with no close neighbors as singleton groups
    connected = set(distance_dict.keys())
    for idx in set(range(n)) - connected:
        groups[idx] = [idx]

    return list(groups.values())


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def combine_verticals(input_folder, output_folder, log_file_path=None, threshold=1.2):
    """Combine vertical GNSS velocities from all .vel files in `input_folder`.

    For each group of stations within `threshold` km of each other, the
    median of the valid (non-NaN) U.vel and U.sig values is written to the
    output file. The first station in the group that has valid U.vel and U.sig
    provides the output Lon, Lat, and Stat.

    Parameters
    ----------
    input_folder : str
        Folder containing 5-column .vel files (Lon Lat U.vel U.sig Stat).
    output_folder : str
        Folder where combined_vertical_velocity_field.vel will be written.
    log_file_path : str or None
        Path for the debug log. Defaults to
        <output_folder>/debug_combined_vertical_velocity_field.log.
    threshold : float
        Proximity threshold in km (default 1.2).

    Returns
    -------
    (output_file_path, log_file_path) : tuple of str
    """
    os.makedirs(output_folder, exist_ok=True)

    if log_file_path is None:
        log_file_path = os.path.join(
            output_folder, 'debug_combined_vertical_velocity_field.log'
        )

    # ------------------------------------------------------------------
    # Load all .vel files
    # ------------------------------------------------------------------
    all_frames = []
    for fname in os.listdir(input_folder):
        if fname.endswith('.vel'):
            fpath = os.path.join(input_folder, fname)
            df = pd.read_csv(
                fpath, sep=r'\s+', header=None,
                names=['Lon', 'Lat', 'U.vel', 'U.sig', 'Stat']
            )
            all_frames.append(df)

    all_data = pd.concat(all_frames, ignore_index=True)
    all_data['U.vel'] = pd.to_numeric(all_data['U.vel'], errors='coerce')
    all_data['U.sig'] = pd.to_numeric(all_data['U.sig'], errors='coerce')

    # ------------------------------------------------------------------
    # Build neighbor groups with BallTree + Union-Find
    # ------------------------------------------------------------------
    stations      = all_data[['Lon', 'Lat']].values
    distance_dict = create_distance_dict(stations, threshold=threshold)
    groups        = make_groups(distance_dict, len(stations))

    # ------------------------------------------------------------------
    # Pre-extract arrays — eliminates per-group pandas overhead
    # ------------------------------------------------------------------
    lon_arr   = all_data['Lon'].to_numpy(dtype=np.float64)
    lat_arr   = all_data['Lat'].to_numpy(dtype=np.float64)
    u_vel_arr = all_data['U.vel'].to_numpy(dtype=np.float64)
    u_sig_arr = all_data['U.sig'].to_numpy(dtype=np.float64)
    stat_arr  = all_data['Stat'].to_numpy()

    # ------------------------------------------------------------------
    # Process each group
    # ------------------------------------------------------------------
    merged_rows = []
    log_lines   = []

    for group in groups:
        g = np.asarray(group, dtype=np.intp)
        u_vel_g = u_vel_arr[g]
        u_sig_g = u_sig_arr[g]

        # Drop NaN independently per column, matching the original notebook's
        # u_vel_valid = group_df['U.vel'].dropna() logic.
        u_vel_valid = u_vel_g[~np.isnan(u_vel_g)]
        u_sig_valid = u_sig_g[~np.isnan(u_sig_g)]

        if len(u_vel_valid) == 0 or len(u_sig_valid) == 0:
            continue  # skip groups with no valid data

        # First valid entry (both U.vel and U.sig non-NaN) provides coordinates
        # and station name, matching original idxmax() on the joint mask.
        both_valid = ~(np.isnan(u_vel_g) | np.isnan(u_sig_g))
        first_idx  = g[both_valid][0] if both_valid.any() else g[0]
        lon          = lon_arr[first_idx]
        lat          = lat_arr[first_idx]
        stat         = stat_arr[first_idx]
        # Apply round() directly to np.float64 scalars to match the original
        # notebook's round(Series.median(), 2) behavior — converting to Python
        # float() first gives a different result at x.xx5 boundaries.
        u_vel_median = round(np.median(u_vel_valid), 2)
        u_sig_median = round(np.median(u_sig_valid), 2)

        merged_rows.append([lon, lat, u_vel_median, u_sig_median, stat])

        # Accumulate log lines — written in one batch at the end
        log_lines.append('Group:')
        for idx in g:
            log_lines.append(
                f"{lon_arr[idx]} {lat_arr[idx]} {u_vel_arr[idx]} {u_sig_arr[idx]} {stat_arr[idx]}"
            )
        log_lines.append(f"Merged: {lon} {lat} {u_vel_median} {u_sig_median} {stat}")
        log_lines.append('')

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    merged_df = pd.DataFrame(merged_rows, columns=['Lon', 'Lat', 'U.vel', 'U.sig', 'Stat'])
    output_file_path = os.path.join(output_folder, 'combined_vertical_velocity_field.vel')
    merged_df.to_csv(output_file_path, sep=' ', index=False, header=False)

    with open(log_file_path, 'w') as f:
        f.write('\n'.join(log_lines))
        f.write('\n')  # trailing newline matches original line-by-line write behavior

    return output_file_path, log_file_path


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python combine_verticals.py <input_folder> <output_folder>')
        sys.exit(1)

    input_folder  = sys.argv[1]
    output_folder = sys.argv[2]

    t0 = time.time()
    out_vel, out_log = combine_verticals(input_folder, output_folder)
    elapsed = time.time() - t0

    print(f"Combined vertical velocity field saved to: {out_vel}")
    print(f"Debug log saved to: {out_log}")
    print(f"Time taken: {elapsed:.2f}s")
