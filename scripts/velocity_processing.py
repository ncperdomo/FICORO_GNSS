"""Scientific velocity-field processing: filtering, uncertainty bounds, manual
outlier removal.

Functions
---------
filter_postseismic_stations      Remove stations near earthquake epicentres.
apply_uncertainty_lower_bound    Enforce a minimum value on velocity uncertainties.
apply_manual_filter              Apply geographic radius-based removal to combined
                                 velocity fields; runs in parallel across all frames.
filter_vertical_velocities_by_sigma  Remove vertical velocities > N σ from mean.
"""

import glob
import os
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from geo_utils import haversine

# Standard column order for 13-column velocity files
_VEL_COLS = [
    'Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj',
    'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat',
]


def filter_postseismic_stations(file_path, log_path, epicentres,
                                 distance_threshold_km):
    """Remove GNSS stations within *distance_threshold_km* of any epicentre.

    Overwrites *file_path* with the filtered data.  Writes removed stations
    to *log_path*.

    Parameters
    ----------
    file_path             : str   Space-separated .vel / .csv file with a header row.
    log_path              : str   Destination for the removed-stations log.
    epicentres            : list  List of (lon, lat) tuples in decimal degrees.
    distance_threshold_km : float Removal radius in kilometres.
    """
    os.makedirs(os.path.dirname(log_path) or '.', exist_ok=True)

    df = pd.read_csv(file_path, sep=r'\s+', header=0)
    original_n = len(df)
    removed = pd.DataFrame(columns=df.columns)

    for lon, lat in epicentres:
        dist = haversine(df['Lon'].values, df['Lat'].values, lon, lat)
        mask = dist <= distance_threshold_km
        removed = pd.concat([removed, df[mask]], ignore_index=True)
        df = df[~mask]

    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.tmp') as tmp:
        df.to_csv(tmp.name, sep=' ', index=False)
        tmp_path = tmp.name
    shutil.move(tmp_path, file_path)
    removed.to_csv(log_path, sep=' ', index=False)

    n_removed = original_n - len(df)
    pct = n_removed / original_n * 100
    print(
        f"Processed {os.path.basename(file_path)}: "
        f"original={original_n}, final={len(df)}, "
        f"removed={n_removed} ({pct:.2f}%)"
    )
    print(f"Log saved to {log_path}")
    print(f"Cleaned data saved to {file_path}")


def apply_uncertainty_lower_bound(csv_folder, file_names, lower_bound=0.5):
    """Enforce *lower_bound* mm/yr on E.sig and N.sig for specific CSV files.

    Files are read from *csv_folder*, modified in place, and saved back.

    Parameters
    ----------
    csv_folder  : str        Folder containing the coherence-filter CSV files.
    file_names  : list[str]  Base filenames (without extension) to process.
    lower_bound : float      Minimum uncertainty value in mm/yr (default 0.5).
    """
    for name in file_names:
        path = os.path.join(csv_folder, name + '.csv')
        if not os.path.exists(path):
            print(f"Warning: {name}.csv not found in {csv_folder}")
            continue

        df = pd.read_csv(path, sep=' ', skiprows=1, header=None)
        df.columns = _VEL_COLS

        df['E.sig'] = df['E.sig'].clip(lower=lower_bound)
        df['N.sig'] = df['N.sig'].clip(lower=lower_bound)

        df.to_csv(path, sep=' ', index=False, header=False)
        print(f"Applied {lower_bound} mm/yr lower bound to uncertainties in {name}.csv")


def apply_manual_filter(temp_combined_folder, combined_folder, filter_criteria_path):
    """Remove GNSS stations inside user-defined geographic exclusion zones.

    Reads every ``combined_vel_*.csv`` from *temp_combined_folder*, applies
    each (center_lon, center_lat, radius_km) criterion from
    *filter_criteria_path*, and writes cleaned files to *combined_folder*.
    Runs all frames in parallel.

    Parameters
    ----------
    temp_combined_folder  : str  Folder with raw combined velocity CSVs.
    combined_folder       : str  Destination for ``*_clean.csv`` and
                                 ``*_removed.log`` files.
    filter_criteria_path  : str  Space-separated CSV with columns:
                                 center_lon, center_lat, radius (km).
    """
    os.makedirs(combined_folder, exist_ok=True)
    filter_criteria = pd.read_csv(filter_criteria_path, sep=' ')
    files = glob.glob(os.path.join(temp_combined_folder, 'combined_vel_*.csv'))

    def _process_one(fpath):
        fname      = os.path.basename(fpath)
        clean_path = os.path.join(combined_folder, fname.replace('.csv', '_clean.csv'))
        log_path   = os.path.join(combined_folder, fname.replace('.csv', '_removed.log'))

        data = pd.read_csv(fpath, sep=' ')
        original_n = len(data)
        log_rows = []

        for _, crit in filter_criteria.iterrows():
            lon_c, lat_c, radius = crit['center_lon'], crit['center_lat'], crit['radius']
            dist = haversine(data['Lon'].values, data['Lat'].values, lon_c, lat_c)
            mask = dist <= radius
            log_rows.append(data[mask])
            data = data[~mask]

            if 'eura' in fname:
                print(f"  ({lon_c}, {lat_c}), r={radius} km → removed {mask.sum()}")

        removed = (pd.concat(log_rows, ignore_index=True) if log_rows
                   else pd.DataFrame(columns=data.columns))
        n_removed = original_n - len(data)

        data.to_csv(clean_path, sep=' ', index=False)
        removed.to_csv(log_path, sep=' ', index=False)

        if 'eura' in fname:
            pct = n_removed / original_n * 100 if original_n else 0
            print(f"  Total removed: {n_removed}/{original_n} ({pct:.2f}%)")

        return clean_path

    with ThreadPoolExecutor() as pool:
        cleaned = list(pool.map(_process_one, files))

    print("--------------------------------------------------------------")
    print("                   All files processed                        ")
    print("--------------------------------------------------------------")
    print("       Check final velocity fields in the paths below:        ")
    print("--------------------------------------------------------------")
    for p in sorted(cleaned):
        print(p)
    print("--------------------------------------------------------------")


def filter_vertical_velocities_by_sigma(input_path, output_path, n_sigma=2):
    """Remove vertical velocities more than *n_sigma* standard deviations from mean.

    Parameters
    ----------
    input_path  : str    Space-separated 5-column file (no header):
                         Lon Lat U.vel U.sig Stat.
    output_path : str    Destination for the filtered file (same format).
    n_sigma     : float  Number of standard deviations for the rejection threshold.

    Returns
    -------
    mean, std : float  Sample mean and std of the *input* U.vel column.
    """
    cols = ['Lon', 'Lat', 'U.vel', 'U.sig', 'Stat']
    df   = pd.read_csv(input_path, sep=r'\s+', header=None, names=cols)

    mean = df['U.vel'].mean()
    std  = df['U.vel'].std()

    keep = (df['U.vel'] >= mean - n_sigma * std) & (df['U.vel'] <= mean + n_sigma * std)
    df[keep].to_csv(output_path, sep=' ', index=False, header=False)
    print(f"Filtered vertical velocities saved to: {output_path}")
    return mean, std
