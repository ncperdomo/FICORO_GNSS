"""File I/O and format-conversion utilities for the FICORO_GNSS pipeline.

Functions
---------
format_raw_velocity_files    Convert .raw files to space-separated .vel files.
convert_csv_to_vel           Convert coherence-filter CSV outputs to .vel format.
prepare_vertical_input_files Extract vertical-velocity columns from full .vel
                             files and stage them in the input_verticals folder.
"""

import os
import shutil

import pandas as pd

# Standard column order for 13-column GNSS velocity files
VEL_COLUMNS = [
    'Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj',
    'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat',
]


def format_raw_velocity_files(input_path, output_path,
                               input_ext='raw', output_ext='vel'):
    """Reformat whitespace-delimited raw velocity files and add a header row.

    Clears *output_path* before writing (creates it if absent).

    Parameters
    ----------
    input_path  : str  Folder containing ``*.{input_ext}`` files.
    output_path : str  Destination folder for ``*.{output_ext}`` files.
    input_ext   : str  Extension of raw input files (default ``'raw'``).
    output_ext  : str  Extension for formatted output files (default ``'vel'``).
    """
    if os.path.exists(output_path):
        for f in os.listdir(output_path):
            fp = os.path.join(output_path, f)
            if os.path.isfile(fp):
                os.remove(fp)
    else:
        os.makedirs(output_path)

    header = ' '.join(VEL_COLUMNS)

    for fname in os.listdir(input_path):
        if not fname.endswith('.' + input_ext):
            continue
        base = os.path.splitext(fname)[0]
        print(f"Translating {base}.{input_ext} into {base}.{output_ext}")
        with open(os.path.join(input_path, fname)) as fh:
            rows = [line.strip().split() for line in fh.readlines()]
        formatted = [header] + [' '.join(row) for row in rows]
        out_path = os.path.join(output_path, f'{base}.{output_ext}')
        with open(out_path, 'w') as fh:
            fh.write('\n'.join(formatted))


def convert_csv_to_vel(csv_folder, vel_folder):
    """Convert space-separated CSV files (from coherence filter) to .vel format.

    Each CSV has a single header row which is stripped.  A leading space is
    prepended to each data line for velrot / pycvframe compatibility.

    Clears *vel_folder* before writing (creates it if absent).

    Parameters
    ----------
    csv_folder : str  Folder containing ``*.csv`` files.
    vel_folder : str  Destination folder for ``*.vel`` files.
    """
    if os.path.exists(vel_folder):
        for f in os.listdir(vel_folder):
            fp = os.path.join(vel_folder, f)
            if os.path.isfile(fp):
                os.remove(fp)
    else:
        os.makedirs(vel_folder)

    for fname in os.listdir(csv_folder):
        if not fname.endswith('.csv'):
            continue
        src = os.path.join(csv_folder, fname)
        dst = os.path.join(vel_folder, os.path.splitext(fname)[0] + '.vel')
        with open(src) as fh:
            lines = [line.strip() for line in fh.readlines()[1:]]  # skip header
        with open(dst, 'w') as fh:
            fh.write('\n'.join(' ' + line for line in lines))


def prepare_vertical_input_files(input_files, destination_folder,
                                  modify_files=None, raw_input_folder=None):
    """Stage vertical-velocity columns from full .vel files for combination.

    For each path in *input_files*:
    - Reads the 13-column .vel file.
    - Keeps only Lon, Lat, U.vel, U.sig, Stat.
    - Replaces 0.00 U.vel / U.sig with NaN for files listed in *modify_files*.
    - Replaces negative U.sig values with their absolute value (warns).
    - Writes a 5-column space-separated file (no header) to *destination_folder*.

    Also copies any ``*.raw`` files from *raw_input_folder* (if it exists) as
    ``*.vel`` files.

    Parameters
    ----------
    input_files        : list of str  Full paths to 13-column .vel sources.
    destination_folder : str          Output folder for 5-column .vel files.
    modify_files       : list of str  Base filenames whose 0.00 values become NaN.
    raw_input_folder   : str or None  Folder with levelling ``*.raw`` files.
    """
    modify_files = modify_files or []
    os.makedirs(destination_folder, exist_ok=True)

    for file_path in input_files:
        if not os.path.exists(file_path):
            print(f"File does not exist: {file_path}")
            continue

        data = pd.read_csv(file_path, sep=r'\s+', header=None, names=VEL_COLUMNS)
        sel = data[['Lon', 'Lat', 'U.vel', 'U.sig', 'Stat']].copy()

        base = os.path.basename(file_path).split('.')[0]
        if base in modify_files:
            sel.loc[sel['U.vel'] == 0.00, 'U.vel'] = pd.NA
            sel.loc[sel['U.sig'] == 0.00, 'U.sig'] = pd.NA

        neg_mask = sel['U.sig'].notna() & (sel['U.sig'] < 0)
        if neg_mask.any():
            print(f"Warning: negative U.sig in {base}.vel — replacing with |value|")
            sel.loc[neg_mask, 'U.sig'] = sel.loc[neg_mask, 'U.sig'].abs()

        out = os.path.join(destination_folder, os.path.basename(file_path))
        sel.to_csv(out, sep=' ', index=False, header=False)

    if raw_input_folder and os.path.exists(raw_input_folder):
        for rf in os.listdir(raw_input_folder):
            if rf.endswith('.raw'):
                src = os.path.join(raw_input_folder, rf)
                dst = os.path.join(destination_folder, rf.replace('.raw', '.vel'))
                shutil.copy(src, dst)
    elif raw_input_folder:
        print(f"Folder does not exist: {raw_input_folder}")
