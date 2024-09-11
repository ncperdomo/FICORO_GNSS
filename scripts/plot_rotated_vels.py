import os
import sys
import glob
import pandas as pd
import numpy as np
import pygmt
from geopy.distance import geodesic

def plot_gps_velocity_fields(folder_path, figure_folder):
    # Find all CSV files in the output_coherence_analysis folder
    file_names = glob.glob(os.path.join(folder_path, '*.csv'))

    # Create a new figure for each file
    for file_name in file_names:
        # Create a new figure
        fig = pygmt.Figure()

        # Set the region and projection of the map
        #fig.basemap(region=[-15, 70, 5, 60], projection='M10c', frame='afg') # Only Euromediterranean and Middle East regions
        fig.basemap(region=[-20, 125, 5, 60], projection='M20c', frame='af') # Entire Alpine-Himalayan belt

        # Create a custom color palette for the relief shading
        pygmt.makecpt(cmap="gray95,gray90,gray85", series=[-10000, 10000, 100])

        # Add shaded topography with transparency
        fig.grdimage(grid="@earth_relief_03m", cmap=True, shading=True, transparency=20) # nan_transparent=True results in error in the latest GMT version

        # Add coastlines
        fig.coast(water='white', borders="1/0.1p,gray90", shorelines="0.1p,black", area_thresh=4000, resolution='h')

        # Read the CSV file from output_coherence_analysis folder
        #df = pd.read_csv(file_name, sep='\s+', skiprows=1, header=None)
        #df.columns = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat']
        usecols = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat'] # I do this to ignore potential extra columns in the CSV file
        df = pd.read_csv(file_name, sep='\s+', skiprows=1, header=None, usecols=range(len(usecols)))
        df.columns = usecols

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

        # Normalise the velocity magnitude to the range [0, 1]
        normalized_vel_mag = (vel_mag - vel_mag.min()) / (vel_mag.max() - vel_mag.min())

        # Create a list to store the vectors
        vectors = []

        # Iterate over each site
        for i in range(len(df)):
            x_start = lon[i]
            y_start = lat[i]
            direction_degrees = np.degrees(np.arctan2(n_vel[i], e_vel[i]))
            length = normalized_vel_mag[i] #* 0.5

            # Add the vector to the list
            vectors.append([x_start, y_start, direction_degrees, length])

        # Plot the GPS velocity vectors from output_coherence_analysis folder (blue)
        fig.plot(
            style='v0.1c+e+n0.15',
            data=vectors,
            fill='red',
            pen='black',
            label='Accepted vel.',
        )

        # Add scale bar
        with pygmt.config(FONT_ANNOT_PRIMARY='8p', FONT_LABEL='8p'):
            fig.basemap(map_scale="JBR+o-9c/-0.8c+c0+w1000k+f+lkm")

        # Add scale vectors to the plot
        # Define origin of the scale vectors
        scale_origin_lon = 68  # Longitude of the scale vector origin
        scale_origin_lat = 16  # Latitude of the scale vector origin

        # Define scale vectors (30 mm/yr normalised appropriately)
        scale_vector_length = 30  # in mm/yr, this will be normalised

        # Calculate the normalised length for the scale vectors
        normalized_scale_length = (scale_vector_length - vel_mag.min()) / (vel_mag.max() - vel_mag.min())

        scale_vectors = [
            [scale_origin_lon, scale_origin_lat, 0, normalized_scale_length],  # Eastward vector
            [scale_origin_lon, scale_origin_lat, 90, normalized_scale_length]  # Northward vector
        ]

        # Plot the scale vectors
        fig.plot(
            style='v0.1c+e+n0.15',
            data=scale_vectors,
            fill='red',
            pen='black',
            label='Accepted vel.',
        )

        # Annotate the scale vectors
        fig.text(
            text=f'{scale_vector_length} mm/yr',
            x=scale_origin_lon - 5,
            y=scale_origin_lat,
            font='7p,black'
        )

        # Get the base file name without extension
        base_name = os.path.splitext(os.path.basename(file_name))[0]

        print(f"Plotting GPS velocities for {base_name} data set")

        # Show the figure
        #fig.show()

        # Create a directory to store the figure files
        os.makedirs(figure_folder, exist_ok=True)

        # Save the figure
        figure_file_pdf = os.path.join(figure_folder, f'{base_name}_map.pdf')
        #figure_file_jpg = os.path.join(figure_folder, f'{base_name}_map.jpg')
        #figure_file_png = os.path.join(figure_folder, f'{base_name}_map.png')

        fig.savefig(figure_file_pdf, dpi=300)
        #fig.savefig(figure_file_jpg, dpi=300)
        #fig.savefig(figure_file_png, dpi=300)

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python plot_rotated_vels.py ./path2/input_folder")
        sys.exit(1)

    input_folder = sys.argv[1]
    figures_path = './results/figures'
    plot_gps_velocity_fields(input_folder, figures_path)
