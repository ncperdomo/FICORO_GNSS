import os
import sys
import glob
import pandas as pd
import numpy as np
import pygmt

_DEFAULT_REGION = [-20, 125, 5, 60]


def plot_gps_velocity_fields(folder_path, figure_folder, region=None):
    """Plot rotated GPS velocity vectors for each CSV in folder_path.

    Parameters
    ----------
    folder_path   : str        Folder of rotated/combined CSV velocity files.
    figure_folder : str        Destination folder for output PDFs.
    region        : list|None  PyGMT region [lon_min, lon_max, lat_min, lat_max].
    """
    _region = region or _DEFAULT_REGION
    file_names = glob.glob(os.path.join(folder_path, '*.csv'))
    os.makedirs(figure_folder, exist_ok=True)

    usecols = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj',
               'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat']

    for file_name in file_names:
        df = pd.read_csv(file_name, sep=r'\s+', skiprows=1, header=None,
                         usecols=range(len(usecols)))
        df.columns = usecols

        if df.shape[0] == 0:
            print(f"Skipping empty file: {file_name}")
            continue

        lon   = df['Lon'].to_numpy()
        lat   = df['Lat'].to_numpy()
        e_vel = df['E.vel'].to_numpy()
        n_vel = df['N.vel'].to_numpy()

        vel_mag = np.sqrt(e_vel**2 + n_vel**2)
        normalized_vel_mag = (vel_mag - vel_mag.min()) / (vel_mag.max() - vel_mag.min())
        direction = np.degrees(np.arctan2(n_vel, e_vel))
        vectors = np.column_stack([lon, lat, direction, normalized_vel_mag]).tolist()

        scale_origin_lon = 68
        scale_origin_lat = 16
        scale_vector_length = 30
        normalized_scale_length = (scale_vector_length - vel_mag.min()) / (vel_mag.max() - vel_mag.min())
        scale_vectors = [
            [scale_origin_lon, scale_origin_lat, 0,  normalized_scale_length],
            [scale_origin_lon, scale_origin_lat, 90, normalized_scale_length],
        ]

        fig = pygmt.Figure()
        fig.basemap(region=_region, projection='M20c', frame='af')
        pygmt.makecpt(cmap="gray95,gray90,gray85", series=[-10000, 10000, 100])
        fig.grdimage(grid="@earth_relief_03m", cmap=True, shading=True, transparency=20)
        fig.coast(water='white', borders="1/0.1p,gray90",
                  shorelines="0.1p,black", area_thresh=4000, resolution='h')

        fig.plot(style='v0.1c+e+n0.15', data=vectors, fill='red', pen='black',
                 label='Accepted vel.')
        fig.plot(style='v0.1c+e+n0.15', data=scale_vectors, fill='red', pen='black',
                 label='Accepted vel.')
        fig.text(text=f'{scale_vector_length} mm/yr',
                 x=scale_origin_lon - 5, y=scale_origin_lat, font='7p,black')

        with pygmt.config(FONT_ANNOT_PRIMARY='8p', FONT_LABEL='8p'):
            fig.basemap(map_scale="JBR+o-9c/-0.8c+c0+w1000k+f+lkm")

        base_name = os.path.splitext(os.path.basename(file_name))[0]
        print(f"Plotting GPS velocities for {base_name}")
        fig.savefig(os.path.join(figure_folder, f'{base_name}_map.pdf'), dpi=300)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_rotated_vels.py ./path/to/input_folder")
        sys.exit(1)

    plot_gps_velocity_fields(sys.argv[1], './results/figures')
