import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import lognorm
import time
import warnings

warnings.simplefilter("ignore", category=RuntimeWarning)
plt.rcParams['figure.max_open_warning'] = 50


def filter_and_plot_data(folder_path, log_output_folder, output_folder, figure_folder):
    """Fit a lognormal distribution to velocity uncertainties and remove
    stations that exceed the 99th percentile.

    Parameters
    ----------
    folder_path       : str  Folder containing ``*.vel`` files.
    log_output_folder : str  Destination for excluded-station CSV files.
    output_folder     : str  Destination for filtered CSV files.
    figure_folder     : str  Destination for diagnostic PDF figures.
    """
    file_names = glob.glob(os.path.join(folder_path, '*.vel'))

    os.makedirs(log_output_folder, exist_ok=True)
    os.makedirs(output_folder,     exist_ok=True)
    os.makedirs(figure_folder,     exist_ok=True)

    for file_name in file_names:
        df = pd.read_csv(file_name, sep=r'\s+', header=0)
        df['E.sig'] = pd.to_numeric(df['E.sig'], errors='coerce')
        df['N.sig'] = pd.to_numeric(df['N.sig'], errors='coerce')

        pos_e = df['E.sig'] > 0
        pos_n = df['N.sig'] > 0

        e_params = lognorm.fit(df.loc[pos_e, 'E.sig'].dropna())
        n_params = lognorm.fit(df.loc[pos_n, 'N.sig'].dropna())

        e_99 = lognorm.ppf(0.99, *e_params)
        n_99 = lognorm.ppf(0.99, *n_params)

        excluded = pd.concat([
            df[df['E.sig'] > e_99],
            df[df['N.sig'] > n_99],
        ]).drop_duplicates()

        filtered_df = df[(df['E.sig'] < e_99) & (df['N.sig'] < n_99)]

        basename  = os.path.splitext(os.path.basename(file_name))[0]
        n_removed = len(excluded)
        n_total   = len(df)
        pct       = n_removed / n_total * 100
        print(f"{basename}: removed {n_removed}/{n_total} ({pct:.2f}%)")

        excluded.to_csv(
            os.path.join(log_output_folder, f'{basename}.csv'), sep=' ', index=False
        )
        filtered_df.to_csv(
            os.path.join(output_folder, f'{basename}.csv'), sep=' ', index=False
        )

        _plot_subfigures(df, basename, figure_folder, e_99, n_99, e_params, n_params)


def _plot_subfigures(df, basename, figure_folder, e_99, n_99, e_params, n_params):
    e_vals = df['E.sig'].dropna()
    n_vals = df['N.sig'].dropna()

    x_e = np.linspace(e_vals.min(), e_vals.max(), 1000)
    x_n = np.linspace(n_vals.min(), n_vals.max(), 1000)
    y_e = lognorm.pdf(x_e, *e_params)
    y_n = lognorm.pdf(x_n, *n_params)

    e_mean = e_vals.mean()
    n_mean = n_vals.mean()

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))

    for ax, vals, x, y, p99, mean, label in [
        (axs[0], e_vals, x_e, y_e, e_99, e_mean, 'East'),
        (axs[1], n_vals, x_n, y_n, n_99, n_mean, 'North'),
    ]:
        counts, bins, _ = ax.hist(vals, bins=20, density=False, alpha=0.7)
        bw = bins[1] - bins[0]
        ax.plot(x, y * len(vals) * bw, 'r-', label='Lognormal fit')
        ax.axvline(p99,  color='g',      linestyle='--', label=f'99 %: {p99:.2f}')
        ax.axvline(mean, color='orange', linestyle='--', label=f'Mean: {mean:.2f}')
        ax.set_title(f'{basename}: {label} velocity uncertainty')
        ax.set_xlabel(f'{label} velocity uncertainty (mm/yr)')
        ax.set_ylabel('Counts')
        ax.set_xlim(-0.5, 6)
        ax.legend()

    plt.tight_layout()
    plt.savefig(
        os.path.join(figure_folder, f'{basename}_lognorm_filter.pdf'),
        dpi=300, format='pdf',
    )
    plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python lognorm_filter.py <input_folder> <output_folder> "
              "<log_output_folder> <figure_folder>")
        sys.exit(1)

    t0 = time.time()
    filter_and_plot_data(sys.argv[1], sys.argv[3], sys.argv[2], sys.argv[4])
    print(f"lognorm_filter done in {time.time() - t0:.2f} s")
