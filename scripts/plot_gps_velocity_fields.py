"""Plot horizontal GPS velocity fields from combined .csv files.

Each CSV is space-separated with 4 header rows followed by data rows:
    Lon Lat E.vel N.vel E.adj N.adj E.sig N.sig Corr U.vel U.adj U.sig Stat

Public API
----------
plot_gps_velocity_fields(folder_path, plate_name=None)
    Plot all velocity fields in folder_path, or just the one whose filename
    contains plate_name (case-insensitive).
"""

import glob
import os

import numpy as np
import pandas as pd
import pygmt


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _select_files(folder_path, plate_name):
    """Return the list of CSV files to plot.

    Parameters
    ----------
    folder_path : str
        Directory that contains ``*.csv`` velocity-field files.
    plate_name : str or None
        If None, return all ``*.csv`` files in the folder.
        Otherwise return only files whose basename contains ``plate_name``
        (case-insensitive).  Raises ``FileNotFoundError`` when a non-None
        plate_name matches zero files.
    """
    all_files = sorted(glob.glob(os.path.join(folder_path, "*.csv")))

    if plate_name is None:
        return all_files

    key = plate_name.lower()
    matches = [f for f in all_files if key in os.path.basename(f).lower()]

    if not matches:
        available = [os.path.splitext(os.path.basename(f))[0] for f in all_files]
        raise FileNotFoundError(
            f"No CSV file in '{folder_path}' contains '{plate_name}'. "
            f"Available files: {available}"
        )

    return matches


def _load_file(file_name):
    """Read a velocity-field CSV.  Returns the DataFrame, or None if empty."""
    df = pd.read_csv(file_name, sep=r"\s+", skiprows=4, header=None)
    df.columns = [
        "Lon", "Lat", "E.vel", "N.vel", "E.adj", "N.adj",
        "E.sig", "N.sig", "Corr", "U.vel", "U.adj", "U.sig", "Stat",
    ]
    if df.shape[0] == 0:
        return None
    return df


def _build_vectors(df):
    """Return (vectors, vel_mag) for the velocity data in df.

    vectors  : list of [lon, lat, direction_deg, normalised_length]
    vel_mag  : 1-D array of raw magnitudes (mm/yr), used for scale vector
    """
    e_vel = df["E.vel"].to_numpy()
    n_vel = df["N.vel"].to_numpy()
    lon   = df["Lon"].to_numpy()
    lat   = df["Lat"].to_numpy()

    vel_mag   = np.sqrt(e_vel**2 + n_vel**2)
    mag_min, mag_max = vel_mag.min(), vel_mag.max()
    norm_mag  = (vel_mag - mag_min) / (mag_max - mag_min)
    direction = np.degrees(np.arctan2(n_vel, e_vel))

    vectors = np.column_stack([lon, lat, direction, norm_mag]).tolist()
    return vectors, vel_mag


def _plot_figure(file_name, df):
    """Create and display one PyGMT figure for the velocity field in df."""
    vectors, vel_mag = _build_vectors(df)

    fig = pygmt.Figure()

    # Base map
    fig.basemap(region=[-20, 125, 5, 60], projection="M20c", frame="af")

    # Shaded topography
    pygmt.makecpt(cmap="gray95,gray90,gray85", series=[-10000, 10000, 100])
    fig.grdimage(grid="@earth_relief_03m", cmap=True, shading=True, transparency=20)

    # Coastlines
    fig.coast(
        water="white",
        borders="1/0.1p,gray90",
        shorelines="0.1p,black",
        area_thresh=4000,
        resolution="h",
    )

    # Velocity vectors
    fig.plot(
        style="v0.1c+e+n0.15",
        data=vectors,
        fill="red",
        pen="black",
        label="Accepted vel.",
    )

    # Scale bar
    with pygmt.config(FONT_ANNOT_PRIMARY="8p", FONT_LABEL="8p"):
        fig.basemap(map_scale="JBR+o-9c/-0.8c+c0+w1000k+f+lkm")

    # Reference scale vectors (30 mm/yr)
    scale_origin_lon = 68
    scale_origin_lat = 16
    scale_vector_length = 30  # mm/yr
    mag_min, mag_max = vel_mag.min(), vel_mag.max()
    norm_scale = (scale_vector_length - mag_min) / (mag_max - mag_min)

    scale_vectors = [
        [scale_origin_lon, scale_origin_lat, 0,  norm_scale],   # Eastward
        [scale_origin_lon, scale_origin_lat, 90, norm_scale],   # Northward
    ]
    fig.plot(
        style="v0.1c+e+n0.15",
        data=scale_vectors,
        fill="red",
        pen="black",
        label="Accepted vel.",
    )
    fig.text(
        text=f"{scale_vector_length} mm/yr",
        x=scale_origin_lon - 5,
        y=scale_origin_lat,
        font="7p,black",
    )

    base_name = os.path.splitext(os.path.basename(file_name))[0]
    print(f"Plotting GPS velocities: {base_name}")
    fig.show()


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def plot_gps_velocity_fields(folder_path, plate_name=None):
    """Plot horizontal GPS velocity fields from CSV files in folder_path.

    Parameters
    ----------
    folder_path : str
        Directory containing space-separated ``*.csv`` velocity-field files
        (13-column format with 4 header rows).
        Default location used in the notebook:
        ``'./results/combined_velocities/manual_filter/'``
    plate_name : str or None, optional
        If None (default), all ``*.csv`` files in folder_path are plotted.
        If a string, only the file(s) whose name contains that string
        (case-insensitive) are plotted.  Examples: ``'eura'``, ``'arab'``,
        ``'igb14'``.
        Raises ``FileNotFoundError`` if no matching file is found.

    Examples
    --------
    # Plot every reference frame:
    plot_gps_velocity_fields('./results/combined_velocities/manual_filter/')

    # Plot only the Eurasia-fixed field:
    plot_gps_velocity_fields(
        './results/combined_velocities/manual_filter/',
        plate_name='eura',
    )
    """
    files = _select_files(folder_path, plate_name)

    for file_name in files:
        df = _load_file(file_name)
        if df is None:
            print(f"Skipping empty file: {file_name}")
            continue
        _plot_figure(file_name, df)
