import argparse
import json
import os
import sys
import glob
import threading
import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree
import concurrent.futures
import time

# ── thread-safe console output ───────────────────────────────────────────────
_print_lock = threading.Lock()

def _tprint(*args, **kwargs):
    """Print with a lock so parallel workers do not interleave lines."""
    with _print_lock:
        print(*args, **kwargs)


# ── per-location sigma level ─────────────────────────────────────────────────
def get_region_stringency(lon, lat, regions):
    """Return the sigma threshold for (lon, lat) given the user's region list."""
    if not regions:
        return 2
    for region in regions:
        if (region['min_lon'] <= lon <= region['max_lon'] and
                region['min_lat'] <= lat <= region['max_lat']):
            return region['sigma']
    return 2


# ── single-file coherence filter ─────────────────────────────────────────────
def filter_gps_velocities(file_name, output_clean_folder, output_excluded_folder,
                           radius=20, geo_strict=False, regions=None,
                           special_case_file=None):
    """Remove velocity outliers using a local Z-score test within *radius* km.

    For each station, the BallTree finds all neighbours within *radius* km in
    O(N log N).  If ≥ 5 neighbours exist, stations whose E.vel or N.vel
    deviates more than *sigma* standard deviations from the local mean are
    flagged as outliers.

    Parameters
    ----------
    file_name              : str   Space-separated CSV with a header row.
    output_clean_folder    : str   Destination for cleaned station CSV.
    output_excluded_folder : str   Destination for excluded station CSV.
    radius                 : float Neighbourhood radius in km (default 20).
    geo_strict             : bool  Apply region-specific sigma levels if True.
    regions                : list  Region dicts with min/max_lon/lat and sigma.
    special_case_file      : str   Basename pattern — skip filtering for matches.
    """
    regions = regions or []
    os.makedirs(output_clean_folder,    exist_ok=True)
    os.makedirs(output_excluded_folder, exist_ok=True)

    df = pd.read_csv(file_name, sep=' ', skiprows=1, header=None)
    df.columns = [
        'Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj',
        'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat',
    ]

    basename          = os.path.splitext(os.path.basename(file_name))[0]
    filtered_stations = set()

    if not (special_case_file and special_case_file in file_name):
        lon_vals = df['Lon'].to_numpy(dtype=float)
        lat_vals = df['Lat'].to_numpy(dtype=float)
        e_vel    = df['E.vel'].to_numpy(dtype=float)
        n_vel    = df['N.vel'].to_numpy(dtype=float)

        # Build BallTree once — O(N log N) neighbourhood search instead of O(N²)
        coords_rad   = np.radians(np.column_stack([lat_vals, lon_vals]))
        radius_rad   = radius / 6371.0
        tree         = BallTree(coords_rad, metric='haversine')
        nbr_idx_list = tree.query_radius(coords_rad, r=radius_rad)

        for i, nbr_idx in enumerate(nbr_idx_list):
            if len(nbr_idx) < 5:
                continue
            sigma = (get_region_stringency(lon_vals[i], lat_vals[i], regions)
                     if geo_strict else 2)

            nbr_e = e_vel[nbr_idx]
            nbr_n = n_vel[nbr_idx]
            # ddof=1 matches pandas Series.std() used in the original implementation
            e_mean, e_std = nbr_e.mean(), nbr_e.std(ddof=1)
            n_mean, n_std = nbr_n.mean(), nbr_n.std(ddof=1)

            outlier_mask = (
                (nbr_e < e_mean - sigma * e_std) |
                (nbr_e > e_mean + sigma * e_std) |
                (nbr_n < n_mean - sigma * n_std) |
                (nbr_n > n_mean + sigma * n_std)
            )
            # Convert positional BallTree indices → DataFrame index labels
            filtered_stations.update(df.index[nbr_idx[outlier_mask]])

    excluded_df = df.loc[sorted(filtered_stations)].drop_duplicates()
    included_df = df.drop(sorted(filtered_stations)).drop_duplicates()

    excluded_df.to_csv(
        os.path.join(output_excluded_folder, f'{basename}.csv'), sep=' ', index=False
    )
    included_df.to_csv(
        os.path.join(output_clean_folder, f'{basename}.csv'), sep=' ', index=False
    )

    n_removed = len(excluded_df)
    n_total   = len(df)
    pct       = n_removed / n_total * 100
    _tprint(f"{basename}: removed {n_removed}/{n_total} ({pct:.2f}%)")


# ── parallel entry point ──────────────────────────────────────────────────────
def parallel_filter_gps_velocities(folder_path, output_clean_folder,
                                    output_excluded_folder, radius=20,
                                    geo_strict=False, regions=None,
                                    special_case_file=None):
    """Run :func:`filter_gps_velocities` in parallel across all CSV files.

    Parameters
    ----------
    folder_path            : str   Folder containing ``*.csv`` files.
    output_clean_folder    : str   Destination for cleaned station CSVs.
    output_excluded_folder : str   Destination for excluded station CSVs.
    radius                 : float Neighbourhood radius in km (default 20).
    geo_strict             : bool  Apply region-specific sigma levels if True.
    regions                : list  Region dicts (see :func:`filter_gps_velocities`).
    special_case_file      : str   Basename pattern to skip filtering for.
    """
    regions    = regions or []
    file_names = glob.glob(os.path.join(folder_path, '*.csv'))

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(
                filter_gps_velocities, fn,
                output_clean_folder, output_excluded_folder,
                radius, geo_strict, regions, special_case_file,
            )
            for fn in file_names
        ]
        concurrent.futures.wait(futures)
        for f in futures:
            f.result()   # re-raise any worker exceptions


# ── CLI ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Coherence filter for GNSS velocity fields.'
    )
    parser.add_argument('folder_path',
                        help='Folder containing lognorm-filtered CSV files')
    parser.add_argument('--output_clean',
                        default='./results/output_coherence_analysis',
                        help='Output folder for cleaned CSVs')
    parser.add_argument('--output_excluded',
                        default='./results/sites_excluded_coherence',
                        help='Output folder for excluded-station CSVs')
    parser.add_argument('--geo_strict',        action='store_true')
    parser.add_argument('--regions_json',      default='')
    parser.add_argument('--special_case_file', default='')
    args = parser.parse_args()

    regions = []
    if args.geo_strict and args.regions_json:
        with open(args.regions_json) as fh:
            regions = json.load(fh)

    t0 = time.time()
    parallel_filter_gps_velocities(
        args.folder_path,
        args.output_clean, args.output_excluded,
        geo_strict=args.geo_strict, regions=regions,
        special_case_file=args.special_case_file,
    )
    print(f"coherence_filter done in {time.time() - t0:.2f} s")
