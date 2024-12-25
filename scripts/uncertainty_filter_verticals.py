import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import lognorm

class UncertaintyFilterVerticals:
    """This class is designed to filter out vertical velocities with uncertainties beyond the 99% of the fitted lognormal distribution."""
    def __init__(self, input_file, output_folder, figures_path):
        self.input_file = input_file
        self.output_folder = output_folder
        self.figures_path = figures_path
        os.makedirs(self.output_folder, exist_ok=True)
        os.makedirs(self.figures_path, exist_ok=True)

    def read_vertical_velocities(self):
        """
        Reads the vertical velocity file and returns a DataFrame.
        """
        col_names = ['Lon', 'Lat', 'U.vel', 'U.sig', 'Stat']
        self.data = pd.read_csv(self.input_file, delim_whitespace=True, header=None, names=col_names)
        return self.data

    def plot_uncertainty_distribution(self, save_fig=True):
        """
        Plots the uncertainty distribution for vertical velocities.
        """
        positive_uncertainties = self.data['U.sig'][self.data['U.sig'] > 0]

        # Fit a lognormal distribution
        shape, loc, scale = lognorm.fit(positive_uncertainties, floc=0)

        # Create a range of x values for plotting the lognormal PDF
        x_vals = np.linspace(positive_uncertainties.min(), positive_uncertainties.max(), 10000)
        pdf_vals = lognorm.pdf(x_vals, shape, loc=loc, scale=scale)

        # Compute the 99th percentile
        p99 = lognorm.ppf(0.99, shape, loc=loc, scale=scale)

        # Plot histogram and fitted distribution
        plt.figure(figsize=(10, 6))
        plt.hist(positive_uncertainties, bins=350, alpha=0.7, color='lightgray', edgecolor='black', density=True, label='Vertical velocity (density)')
        plt.plot(x_vals, pdf_vals, label='Lognormal fit', color='red')
        plt.axvline(p99, color='black', linestyle='dashed', linewidth=1.5, label=f'99th percentile: {p99:.2f}')
        plt.xlabel('Vertical velocity uncertainty (mm/yr)')
        plt.ylabel('Density')
        # Crop x-axis
        plt.xlim(0, 6)
        plt.legend()
        plt.tight_layout()

        if save_fig:
            output_figure_path = os.path.join(self.figures_path, 'vertical_uncertainty_distribution.pdf')
            plt.savefig(output_figure_path, format='pdf', dpi=300)
        plt.show()

    def filter_uncertainties(self):
        """
        Filters out velocities with uncertainties beyond the 99% of the fitted lognormal distribution.
        """
        positive_uncertainties = self.data['U.sig'][self.data['U.sig'] > 0]

        # Fit a lognormal distribution
        shape, loc, scale = lognorm.fit(positive_uncertainties, floc=0)

        # Compute the 99th percentile
        p99 = lognorm.ppf(0.99, shape, loc=loc, scale=scale)

        # Filter data
        self.filtered_data = self.data[self.data['U.sig'] <= p99]

        filtered_file_path = os.path.join(self.output_folder, 'filtered_vertical_velocity_field.vel')
        self.filtered_data.to_csv(filtered_file_path, sep=' ', index=False, header=False)

        print(f"Filtered vertical velocities saved to: {filtered_file_path}")
        return self.filtered_data

# Call the class to filter the vertical velocities with the default parameters and input file
if __name__ == "__main__":
    input_file = './results/combined_velocities/manual_filter/combined_vertical_velocity_field.vel'
    output_folder = './results/combined_velocities/manual_filter/'
    figures_path = './results/figures/'

    filter_verticals = UncertaintyFilterVerticals(input_file, output_folder, figures_path)
    filter_verticals.read_vertical_velocities()
    filter_verticals.plot_uncertainty_distribution()
    filter_verticals.filter_uncertainties()
