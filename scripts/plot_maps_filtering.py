import os
import glob
import pandas as pd
import numpy as np
import pygmt
from geopy.distance import geodesic

def plot_gps_velocities(folder_path, excluded_lognorm, excluded_coherence, figure_folder):
    # Find all CSV files in the output_coherence_analysis folder
    file_names = glob.glob(os.path.join(folder_path, '*.csv'))

    # Create a new figure for each file
    for file_name in file_names:
        # Create a new figure
        fig = pygmt.Figure()

        # Set the region and projection of the map
        #fig.basemap(region=[-15, 70, 5, 60], projection='M10c', frame='afg')
        fig.basemap(region=[-20, 125, 5, 60], projection='M15c', frame='afg') # All Alpine-Himalayan belt

        # Create a custom color palette for the relief shading
        pygmt.makecpt(cmap="gray95,gray90,gray85", series=[-10000, 10000, 100])
        
        # Add shaded topography with transparency
        fig.grdimage(grid="@earth_relief_03m", cmap=True, shading=True, transparency=90, nan_transparent=True)

        # Add coastlines
        fig.coast(water='white', borders="1/0.1p,gray90", shorelines="0.1p,black", area_thresh=4000, resolution='h')

        # Read the CSV file from output_coherence_analysis folder
        df = pd.read_csv(file_name, sep=' ', skiprows=1, header=None)
        df.columns = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat']

        # Check if the data frame is empty
        if df.shape[0] == 0:
            print(f"Skipping empty file: {file_name}")
            continue
        
        # Extract the coordinates, E and N velocity components, and E and N sig from the data frame
        lon = df['Lon']
        lat = df['Lat']
        e_vel = df['E.vel']
        n_vel = df['N.vel']
        e_sig = df['E.sig']
        n_sig = df['N.sig']

        # Calculate the velocity magnitude for scaling
        vel_mag = np.sqrt(e_vel**2 + n_vel**2)

        # Normalize the velocity magnitude to the range [0, 1]
        normalized_vel_mag = (vel_mag - vel_mag.min()) / (vel_mag.max() - vel_mag.min())

        # Create a list to store the vectors
        vectors = []

        # Iterate over each site
        for i in range(len(df)):
            x_start = lon[i]
            y_start = lat[i]
            direction_degrees = 90 - np.arctan2(n_vel[i], e_vel[i]) * 180 / np.pi
            length = normalized_vel_mag[i] * 0.5

            # Add the vector to the list
            vectors.append([x_start, y_start, direction_degrees, length])

        # Plot the GPS velocity vectors from output_coherence_analysis folder (blue)
        fig.plot(
            style='v0.1c+e',
            data=vectors,
            fill='blue',
            pen='black',
            label='Accepted vel.',
        )

        # Get the base file name without extension
        base_name = os.path.splitext(os.path.basename(file_name))[0]

        # Read the CSV file from excluded_lognorm folder
        lognorm_file = os.path.join(excluded_lognorm, f"{base_name}.csv")
        if os.path.exists(lognorm_file):
            try:
                df_lognorm = pd.read_csv(lognorm_file, sep=' ', skiprows=1, header=None, on_bad_lines='skip')
                if df_lognorm.shape[0] == 0:
                    pass
                    #print(f"Skipping empty file: {lognorm_file}")
                else:
                    lon_lognorm = df_lognorm[0]
                    lat_lognorm = df_lognorm[1]
                    e_vel_lognorm = df_lognorm[2]
                    n_vel_lognorm = df_lognorm[3]

                    vectors_lognorm = []
                    for j in range(len(df_lognorm)):
                        x_start_lognorm = lon_lognorm[j]
                        y_start_lognorm = lat_lognorm[j]
                        direction_degrees_lognorm = 90 - np.arctan2(n_vel_lognorm[j], e_vel_lognorm[j]) * 180 / np.pi
                        length_lognorm = normalized_vel_mag[j] * 0.5

                        vectors_lognorm.append([x_start_lognorm, y_start_lognorm, direction_degrees_lognorm, length_lognorm])

                    # Plot the GPS velocity vectors from excluded_lognorm folder (orange)
                    fig.plot(
                        style='v0.1c+e',
                        data=vectors_lognorm,
                        fill='orange',
                        pen='orange',
                        label='Filtered lognorm',
                    )

            except pd.errors.EmptyDataError:
                pass
                #print(f"Skipping empty file: {lognorm_file}")
                  
        # Read the CSV file from excluded_coherence folder
        coherence_file = os.path.join(excluded_coherence, f"{base_name}.csv")
        if os.path.exists(coherence_file):
            try:
                df_coherence = pd.read_csv(coherence_file, sep=' ', skiprows=1, header=None, on_bad_lines='skip')
                if df_coherence.shape[0] == 0:
                    pass
                    #print(f"Skipping empty file: {coherence_file}")
                else:
                    lon_coherence = df_coherence[0]
                    lat_coherence = df_coherence[1]
                    e_vel_coherence = df_coherence[2]
                    n_vel_coherence = df_coherence[3]

                    vectors_coherence = []
                    for k in range(len(df_coherence)):
                        x_start_coherence = lon_coherence[k]
                        y_start_coherence = lat_coherence[k]
                        direction_degrees_coherence = 90 - np.arctan2(n_vel_coherence[k], e_vel_coherence[k]) * 180 / np.pi
                        length_coherence = normalized_vel_mag[k] * 0.5

                        vectors_coherence.append([x_start_coherence, y_start_coherence, direction_degrees_coherence, length_coherence])

                    # Plot the GPS velocity vectors from excluded_coherence folder (red)
                    fig.plot(
                        style='v0.1c+e',
                        data=vectors_coherence,
                        fill='red',
                        pen='red',
                        label='Filtered coherence',
                    )

            except pd.errors.EmptyDataError:
                pass
                #print(f"Skipping empty file: {coherence_file}")

        # Add a legend
        fig.legend(position='JTR+o0.15c/-1.25c', box=True)
        
        # Add scale bar
        with pygmt.config(FONT_ANNOT_PRIMARY='8p', FONT_LABEL='8p'):
            fig.basemap(map_scale="JBR+o-1.6c/-0.8c+c0+w1000k+f+lkm")

        print(f"Plotting GPS velocities for {base_name} data set")

        # Show the figure
        # fig.show()

        # Create a directory to store the figure files
        os.makedirs(figure_folder, exist_ok=True)

        # Save the figure
        figure_file_pdf = os.path.join(figure_folder, f'{base_name}_map.pdf')
        #figure_file_jpg = os.path.join(figure_folder, f'{base_name}_map.jpg')
        #figure_file_png = os.path.join(figure_folder, f'{base_name}_map.png')

        fig.savefig(figure_file_pdf, dpi=300)
        #fig.savefig(figure_file_jpg, dpi=300)
        #fig.savefig(figure_file_png, dpi=300)

# Folder paths containing space-separated CSV files
figures_path = './results/figures'
folder_path = './results/output_coherence_analysis'
excluded_lognorm = './results/sites_excluded_lognorm_99'
excluded_coherence = './results/sites_excluded_coherence'

plot_gps_velocities(folder_path, excluded_lognorm, excluded_coherence, figures_path)
