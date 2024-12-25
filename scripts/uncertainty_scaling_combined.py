# uncertainty_scaling.py

import numpy as np
import scipy.stats as stats
from scipy.stats import lognorm
import os
import pandas as pd
import matplotlib.pyplot as plt

def plot_uncertainty_distributions(original_uncertainties, scaled_uncertainties, component, solution_name):
    """
    Plot the original and scaled uncertainty distributions, lognormal fit, 99th percentile, and mean lines.
    """
    # Ensure we only fit the lognormal on positive values
    positive_uncertainties = original_uncertainties[original_uncertainties > 0]
    positive_uncertainties_scaled = scaled_uncertainties[scaled_uncertainties > 0]

    # Fit a lognormal distribution to the positive uncertainties
    shape, loc, scale = lognorm.fit(positive_uncertainties, floc=0)
    shape_scaled, loc_scaled, scale_scaled = lognorm.fit(positive_uncertainties_scaled, floc=0)
    
    # Create a range of x values for plotting the lognormal PDF
    x_vals = np.linspace(min(positive_uncertainties), max(positive_uncertainties), 1000)
    pdf_vals = lognorm.pdf(x_vals, shape, loc=loc, scale=scale)

    x_vals_scaled = np.linspace(min(positive_uncertainties_scaled), max(positive_uncertainties_scaled), 1000)
    pdf_vals_scaled = lognorm.pdf(x_vals_scaled, shape_scaled, loc=loc_scaled, scale=scale_scaled)
    
    # Compute the 99th percentile for the positive uncertainties
    p99 = lognorm.ppf(0.99, shape, loc=loc, scale=scale)
    p99_scaled = lognorm.ppf(0.99, shape_scaled, loc=loc_scaled, scale=scale_scaled)

    # Calculate the means of the raw and scaled uncertainties
    raw_mean = np.mean(original_uncertainties)
    scaled_mean = np.mean(scaled_uncertainties)

    # Plot histograms for raw and scaled uncertainties (count)
    plt.figure(figsize=(10, 6))
    # Dynamically change the label based on the component
    if component == 'E.sig':
        component_string = '(East)'
        plt.hist(original_uncertainties, bins=30, alpha=0.5, label=f'Original uncertainties {component_string}', density=False)
        plt.hist(scaled_uncertainties, bins=30, alpha=0.5, label=f'Scaled uncertainties {component_string}', density=False)
    elif component == 'N.sig':
        component_string = '(North)'
        plt.hist(original_uncertainties, bins=30, alpha=0.5, label=f'Original uncertainties {component_string}', density=False)
        plt.hist(scaled_uncertainties, bins=30, alpha=0.5, label=f'Scaled uncertainties {component_string}', density=False)

    # Plot the fitted lognormal curve (scaled to match the counts)
    count, bins, _ = plt.hist(positive_uncertainties, bins=30, alpha=0.0)  # Get the bin heights for raw data
    scale_factor = max(count) / max(pdf_vals)  # Scale factor to align the PDF to the counts
    plt.plot(x_vals, pdf_vals * scale_factor, label='Lognormal fit (original)', color='red', linewidth=2)

    count_scaled, bins_scaled, _ = plt.hist(positive_uncertainties_scaled, bins=30, alpha=0.0)  # Get the bin heights for raw data
    scale_factor_scaled = max(count_scaled) / max(pdf_vals_scaled)  # Scale factor to align the PDF to the counts
    plt.plot(x_vals_scaled, pdf_vals_scaled * scale_factor_scaled, label='Lognormal fit (scaled)', color='green', linewidth=2)
    
    # Plot the 99th percentile vertical dashed line for the original uncertainties
    plt.axvline(p99, color='gray', linestyle='dashed', linewidth=2, label=f'99th percentile (original): {p99:.2f}')
    plt.axvline(p99_scaled, color='black', linestyle='dashed', linewidth=2, label=f'99th percentile (scaled): {p99_scaled:.2f}')
    
    # Plot vertical dashed lines for the mean of the raw and scaled uncertainties
    plt.axvline(raw_mean, color='blue', linestyle='dashed', linewidth=2, label=f'Mean (original): {raw_mean:.2f}')
    plt.axvline(scaled_mean, color='orange', linestyle='dashed', linewidth=2, label=f'Mean (scaled): {scaled_mean:.2f}')
    
    # Set x-axis limit to a maximum of 5.5 and leave some space for the ligure label on the top left (for the manuscript)
    plt.xlim(-0.2, 5.5)
    plt.ylim(0, 4200)
    
    # Dynamically change the xlabel based on the component
    if component == 'E.sig':
        plt.xlabel('East velocity uncertainty (mm/yr)')
        filename = "uncertainty_scaling_east_component.pdf"
    elif component == 'N.sig':
        plt.xlabel('North velocity uncertainty (mm/yr)')
        filename = "uncertainty_scaling_north_component.pdf"
    plt.ylabel('Number of GNSS stations')
    plt.legend()
    # Save the figures to the results/figures folder
    plt.savefig(f'./results/figures/{filename}', format='pdf', dpi=300)
    plt.show()

def read_velocity_solution(filename):
    """
    Reads a velocity solution file and returns a DataFrame.
    """
    columns = ['Lon', 'Lat', 'E.vel', 'N.vel', 'E.adj', 'N.adj', 'E.sig', 'N.sig', 'Corr', 'U.vel', 'U.adj', 'U.sig', 'Stat']
    return pd.read_csv(filename, names=columns, sep=' ', skiprows=1, header=None, on_bad_lines='skip')

def remove_nans_and_fit_lognormal(uncertainties):
    """
    Removes NaNs and fits a lognormal distribution to positive values.
    """
    # Remove NaN values
    uncertainties = uncertainties.replace([np.inf, -np.inf], np.nan)
    uncertainties = uncertainties.dropna()  # Remove NaNs
    
    # Fit the lognormal to positive values
    positive_uncertainties = uncertainties[uncertainties > 0]
    return positive_uncertainties

def log_normal_params(data):
    """
    Calculate the mean and std of log-transformed data for log-normal distribution, using only positive values.
    """
    # Remove NaNs and get positive uncertainties
    positive_uncertainties = remove_nans_and_fit_lognormal(data)
    
    log_data = np.log(positive_uncertainties)
    mean = np.mean(log_data)
    std = np.std(log_data)
    return mean, std

def scale_uncertainty(original_uncertainty, original_mean, original_std, reference_mean, reference_std):
    """
    Scale an uncertainty from the original distribution to the reference distribution.
    """
    original_uncertainty = pd.Series([original_uncertainty])
    positive_uncertainty = original_uncertainty[original_uncertainty > 0].values[0]  # Ensure positive values
    
    log_original_uncertainty = np.log(positive_uncertainty)
    percentile = stats.norm.cdf((log_original_uncertainty - original_mean) / original_std)
    scaled_uncertainty = np.exp(reference_mean + reference_std * stats.norm.ppf(percentile))
    return scaled_uncertainty

def remove_outliers_lognormal(df, component):
    """
    Removes velocities outside the 99% of the fitted lognormal distribution for a given component.
    """
    # Remove NaNs and get positive uncertainties
    positive_uncertainties = remove_nans_and_fit_lognormal(df[component])

    # Fit the lognormal distribution and compute the 99th percentile
    component_params = lognorm.fit(positive_uncertainties)
    component_99th = lognorm.ppf(0.99, *component_params)

    # Filter out data points that exceed the 99th percentile
    filtered_df = df[df[component] < component_99th]
    return filtered_df, component_99th

def harmonise_uncertainties(input_folder, reference_filename, output_folder):
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    reference_df = read_velocity_solution(reference_filename)

    # Compute reference (target) distribution parameters
    ref_means, ref_stds = {}, {}
    for component in ['E.sig', 'N.sig']:
        ref_means[component], ref_stds[component] = log_normal_params(reference_df[component])

    # Process each solution file in the input folder
    for solution_file in os.listdir(input_folder):
        if solution_file.endswith('.csv'):
            solution_df = read_velocity_solution(os.path.join(input_folder, solution_file))
            solution_name = solution_file.split('.')[0]
            print(f"Processing {solution_name}")

            for component in ['E.sig', 'N.sig']:
                # Store the raw input data from the CSV file before any preprocessing
                raw_uncertainties = solution_df[component].copy()

                # Remove NaNs and fit to positive values
                processed_uncertainties = remove_nans_and_fit_lognormal(solution_df[component])

                # Calculate the lognormal params from positive values
                original_mean, original_std = log_normal_params(processed_uncertainties)
                
                # Scale uncertainties
                scaled_uncertainties = processed_uncertainties.apply(
                    scale_uncertainty,
                    args=(original_mean, original_std, ref_means[component], ref_stds[component])
                )

                # Here I'll just plot the uncertainty distributions for the solution in Eurasia-fixed reference frame, for the manuscript's supplementary material
                if "eura" in solution_name:                       
                    # Plotting uncertainty distributions (original raw vs scaled) with lognormal fit, mean lines, and 99th percentile
                    plot_uncertainty_distributions(raw_uncertainties, scaled_uncertainties, component, solution_name)
                
                # Store the scaled uncertainties in the solution DataFrame
                solution_df[f'{component}.scaled'] = scaled_uncertainties.round(2)

                # Remove outliers based on the scaled uncertainties
                solution_df, _ = remove_outliers_lognormal(solution_df, f'{component}.scaled')

            # Set U.vel and U.adj to 0.00 before saving the scaled CSV file, as combined vertical velocities are computed separately, as explained in the manuscript 
            solution_df['U.vel'] = 0.00
            solution_df['U.adj'] = 0.00
            solution_df['U.sig'] = 0.00

            # Save the scaled solution to the output folder
            output_file_path = os.path.join(output_folder, f'{solution_name}_scaled.csv')
            solution_df.to_csv(output_file_path, index=False, sep='\t')
