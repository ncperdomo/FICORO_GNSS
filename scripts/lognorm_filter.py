import os
import sys 
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import lognorm
import re
import time
import warnings

# Suppress RuntimeWarnings
warnings.simplefilter("ignore", category=RuntimeWarning)
plt.rcParams['figure.max_open_warning'] = 50  # Avoid warnings

print(f"########## Removing outliers based on lognorm uncertainty distribution ###########")

def filter_and_plot_data(folder_path, log_output_folder, output_folder, figure_folder):
    # Find all .vel files in the folder
    file_names = glob.glob(os.path.join(folder_path, '*.vel'))

    # Empty list to store data frames
    dfs = []

    # Load each .vel file as a data frame
    for file_name in file_names:
        with open(file_name, 'r') as f:
            lines = f.readlines()

        data = []
        for line in lines:
            line_data = re.split(r'\s+', line.strip())
            data.append(line_data)

        # Assign column names and remove first row (text) from each data frame
        df = pd.DataFrame(data)
        df.columns = df.iloc[0]
        df = df.iloc[1:]

        # Convert non-numeric values to NaN
        df['E.sig'] = pd.to_numeric(df['E.sig'], errors='coerce')
        df['N.sig'] = pd.to_numeric(df['N.sig'], errors='coerce')

        dfs.append(df)

    # Create a directory to store the CSV files listing excluded sites
    os.makedirs(log_output_folder, exist_ok=True)

    # Create a directory to store the CSV files listing filtered output
    os.makedirs(output_folder, exist_ok=True)

    # Iterate over each data frame
    for i, df in enumerate(dfs):
        # Filter stations with uncertainties higher than the 99th percentile
        e_sig_higher_than_99 = df.loc[df['E.sig'] > np.percentile(df['E.sig'], 99)]
        n_sig_higher_than_99 = df.loc[df['N.sig'] > np.percentile(df['N.sig'], 99)]

        # Combine the stations for both components
        combined_stations_higher_than_99 = pd.concat([e_sig_higher_than_99, n_sig_higher_than_99]).drop_duplicates()

        print(f"----------------------------------------------------------------------------------")
        
        # Print the number of removed stations for the current dataset
        file_name = os.path.splitext(os.path.basename(file_names[i]))[0]
        num_removed = len(combined_stations_higher_than_99)
        num_total = len(df)
        percentage_removed = (num_removed / num_total) * 100
        print(f"Number of stations removed for {file_name}: {num_removed} / {num_total} ({percentage_removed:.2f}%)")

        # Save the stations with uncertainties larger than 99 percentile to a CSV file
        log_output_file = os.path.join(log_output_folder, f'{file_name}.csv')
        combined_stations_higher_than_99.to_csv(log_output_file, sep=' ', index=False)
        print(f"Sites excluded: {log_output_file}")

        # Filter stations with uncertainties smaller than the 99th percentile
        e_sig_smaller_than_99 = df.loc[df['E.sig'] < np.percentile(df['E.sig'], 99)]
        n_sig_smaller_than_99 = df.loc[df['N.sig'] < np.percentile(df['N.sig'], 99)]

        # Combine the stations for both components
        combined_stations_smaller_than_99 = pd.concat([e_sig_smaller_than_99, n_sig_smaller_than_99]).drop_duplicates()

        # Filter stations whose uncertainties in both East and North components are smaller than the 99th percentile
        combined_stations_smaller_than_99 = combined_stations_smaller_than_99[
            (combined_stations_smaller_than_99['E.sig'] < np.percentile(df['E.sig'], 99)) &
            (combined_stations_smaller_than_99['N.sig'] < np.percentile(df['N.sig'], 99))
        ]

        # Save the stations with uncertainties smaller than 99 percentile to a CSV file
        output_file = os.path.join(output_folder, f'{file_name}.csv')
        combined_stations_smaller_than_99.to_csv(output_file, sep=' ', index=False)
        print(f"Filtered velocities: {output_file}")

        # Plot individual subfigures for each dataset
        plot_subfigures(df, file_name, figure_folder)

def plot_subfigures(df, file_name, figure_folder):
    # Remove NaN values
    e_sig_values = df['E.sig'].dropna()
    n_sig_values = df['N.sig'].dropna()

    # Calculate lognormal parameters for E.sig column
    e_sig_params = lognorm.fit(e_sig_values)
    # Calculate lognormal parameters for N.sig column
    n_sig_params = lognorm.fit(n_sig_values)

    # Generate data points for best lognormal fit
    x_e_sig = np.linspace(e_sig_values.min(), e_sig_values.max(), 100)
    y_e_sig = lognorm.pdf(x_e_sig, *e_sig_params)

    x_n_sig = np.linspace(n_sig_values.min(), n_sig_values.max(), 100)
    y_n_sig = lognorm.pdf(x_n_sig, *n_sig_params)

    # Calculate the 99th percentile of the lognormal distribution
    e_sig_99 = lognorm.ppf(0.99, *e_sig_params)
    n_sig_99 = lognorm.ppf(0.99, *n_sig_params)

    # Calculate the median of the E.sig values
    e_sig_median = np.median(e_sig_values)

    # Calculate the median of the N.sig values
    n_sig_median = np.median(n_sig_values)

    # Create two subplots in one figure
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))

    # Plot histogram and lognormal fit in the first subplot
    counts_e_sig, bins_e_sig, _ = axs[0].hist(e_sig_values, bins=20, density=False, alpha=0.7)
    bin_width_e_sig = bins_e_sig[1] - bins_e_sig[0]
    normalized_y_e_sig = y_e_sig * len(e_sig_values) * bin_width_e_sig
    axs[0].plot(x_e_sig, normalized_y_e_sig, 'r-', label='Lognormal Fit')
    axs[0].axvline(e_sig_99, color='k', linestyle='--', label=f'99%: {e_sig_99:.2f}')
    axs[0].axvline(e_sig_median, color='orange', linestyle='--', label=f'Median: {e_sig_median:.2f}')
    axs[0].set_title(f'{file_name}: East Velocity Uncertainty')
    axs[0].set_xlabel('Uncertainty')
    axs[0].set_ylabel('Count')
    axs[0].legend()

    # Plot histogram and lognormal fit in the second subplot
    counts_n_sig, bins_n_sig, _ = axs[1].hist(n_sig_values, bins=20, density=False, alpha=0.7)
    bin_width_n_sig = bins_n_sig[1] - bins_n_sig[0]
    normalized_y_n_sig = y_n_sig * len(n_sig_values) * bin_width_n_sig
    axs[1].plot(x_n_sig, normalized_y_n_sig, 'r-', label='Lognormal Fit')
    axs[1].axvline(n_sig_99, color='k', linestyle='--', label=f'99%: {n_sig_99:.2f}')
    axs[1].axvline(n_sig_median, color='orange', linestyle='--', label=f'Median: {n_sig_median:.2f}')
    axs[1].set_title(f'{file_name}: North Velocity Uncertainty')
    axs[1].set_xlabel('Uncertainty')
    axs[1].set_ylabel('Count')
    axs[1].legend()

    plt.tight_layout()

    # Create a directory to store the figure files
    os.makedirs(figure_folder, exist_ok=True)

    print(f"Saving figure for {file_name}...")

    # Save the figure in high definition as PDF, JPG, and PNG files
    figure_file_pdf = os.path.join(figure_folder, f'{file_name}_lognorm_filter.pdf')

    plt.savefig(figure_file_pdf, dpi=300, format='pdf')

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 5:
        print("Usage: python lognorm_filter.py ./path2/input_folder ./path2/output_folder ./path2/log_output_folder ./path2/figure_folder")
        sys.exit(1)

    folder_path = sys.argv[1]
    output_folder = sys.argv[2]
    log_output_folder = sys.argv[3]
    figure_folder = sys.argv[4]

    # Time the execution of the combine_velocities function
    start_time = time.time()
    filter_and_plot_data(folder_path, log_output_folder, output_folder, figure_folder)
    end_time = time.time()

    # Calculate and print the elapsed time
    elapsed_time = end_time - start_time
    print("Time taken to run lognorm_filter.py: {:.2f} seconds".format(elapsed_time))
