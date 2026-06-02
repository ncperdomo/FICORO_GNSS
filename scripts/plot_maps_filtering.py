import os
import glob
import pandas as pd
import numpy as np
import pygmt

_DEFAULT_REGION = [-20, 125, 5, 60]


def plot_gps_velocities(folder_path, excluded_lognorm, excluded_coherence,
                        figure_folder, region=None):
    """Plot accepted and filtered GPS velocity vectors for each input file.

    Parameters
    ----------
    folder_path        : str        Folder of coherence-filtered CSV files.
    excluded_lognorm   : str        Folder of lognorm-excluded station CSVs.
    excluded_coherence : str        Folder of coherence-excluded station CSVs.
    figure_folder      : str        Destination folder for output PDFs.
    region             : list|None  PyGMT region [lon_min, lon_max, lat_min, lat_max].
    """
    _region = region or _DEFAULT_REGION
    file_names = glob.glob(os.path.join(folder_path, '*.csv'))
    os.makedirs(figure_folder, exist_ok=True)

    for file_name in file_names:
        fig = pygmt.Figure()
        fig.basemap(region=_region, projection='M20c', frame='af')
        pygmt.makecpt(cmap="gray95,gray90,gray85", series=[-10000, 10000, 100])
        fig.grdimage(grid="@earth_relief_03m", cmap=True, shading=True, transparency=20)
        fig.coast(water='white', borders="1/0.1p,gray90",
                  shorelines="0.1p,black", area_thresh=4000, resolution='h')

        df = pd.read_csv(file_name, sep=' ', skiprows=1, header=None)
        df.columns = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj',
                      'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat']

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

        fig.plot(style='v0.1c+e', data=vectors, fill='blue', pen='black',
                 label='Accepted vel.')

        base_name = os.path.splitext(os.path.basename(file_name))[0]

        # Lognorm-excluded stations
        lognorm_file = os.path.join(excluded_lognorm, f"{base_name}.csv")
        if os.path.exists(lognorm_file):
            try:
                df_ln = pd.read_csv(lognorm_file, sep=' ', skiprows=1,
                                    header=None, on_bad_lines='skip')
                if df_ln.shape[0] > 0:
                    e_ln = df_ln.iloc[:, 2].to_numpy()
                    n_ln = df_ln.iloc[:, 3].to_numpy()
                    vm_ln = np.sqrt(e_ln**2 + n_ln**2)
                    # scale within the accepted-station range, halved for visual distinction
                    norm_ln = (vm_ln - vel_mag.min()) / (vel_mag.max() - vel_mag.min()) * 0.5
                    dir_ln  = np.degrees(np.arctan2(n_ln, e_ln))
                    vecs_ln = np.column_stack([
                        df_ln.iloc[:, 0].to_numpy(),
                        df_ln.iloc[:, 1].to_numpy(),
                        dir_ln, norm_ln,
                    ]).tolist()
                    fig.plot(style='v0.1c+e', data=vecs_ln, fill='orange',
                             pen='orange', label='Filtered lognorm')
            except pd.errors.EmptyDataError:
                pass

        # Coherence-excluded stations
        coherence_file = os.path.join(excluded_coherence, f"{base_name}.csv")
        if os.path.exists(coherence_file):
            try:
                df_co = pd.read_csv(coherence_file, sep=' ', skiprows=1,
                                    header=None, on_bad_lines='skip')
                if df_co.shape[0] > 0:
                    e_co = df_co.iloc[:, 2].to_numpy()
                    n_co = df_co.iloc[:, 3].to_numpy()
                    vm_co = np.sqrt(e_co**2 + n_co**2)
                    norm_co = (vm_co - vel_mag.min()) / (vel_mag.max() - vel_mag.min()) * 0.5
                    dir_co  = np.degrees(np.arctan2(n_co, e_co))
                    vecs_co = np.column_stack([
                        df_co.iloc[:, 0].to_numpy(),
                        df_co.iloc[:, 1].to_numpy(),
                        dir_co, norm_co,
                    ]).tolist()
                    fig.plot(style='v0.1c+e', data=vecs_co, fill='red',
                             pen='red', label='Filtered coherence')
            except pd.errors.EmptyDataError:
                pass

        fig.legend(position='JTR+o0.15c/-1.25c', box=True)
        with pygmt.config(FONT_ANNOT_PRIMARY='8p', FONT_LABEL='8p'):
            fig.basemap(map_scale="JBR+o-9c/-0.8c+c0+w1000k+f+lkm")

        print(f"Plotting GPS velocities for {base_name}")
        fig.savefig(os.path.join(figure_folder, f'{base_name}_map.pdf'), dpi=300)


if __name__ == "__main__":
    plot_gps_velocities(
        folder_path='./results/output_coherence_analysis',
        excluded_lognorm='./results/sites_excluded_lognorm_99',
        excluded_coherence='./results/sites_excluded_coherence',
        figure_folder='./results/figures',
    )
