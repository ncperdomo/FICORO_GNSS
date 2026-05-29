"""Euler-rotation and velocity-combination runners.

Functions
---------
get_available_memory_gb   Return available system RAM in GB.
run_euler_rotations       Rotate all IGB14 files to each plate reference frame.
run_velocity_combination  Combine rotated velocity fields (auto / parallel / sequential).
"""

import concurrent.futures
import os
import shutil
import subprocess
import sys

import psutil


def get_available_memory_gb():
    """Return available system memory in GB."""
    return psutil.virtual_memory().available / (1024 ** 3)


def run_euler_rotations(euler_poles, igb14_folder, igb_nocomb_path,
                        results_path, pycvframe_path):
    """Apply Euler-pole rotations to every IGB14 .vel file.

    Creates one destination folder per plate under
    ``<results_path>/<igb_nocomb_path>/<plate>`` and moves rotated files there.

    Parameters
    ----------
    euler_poles    : dict  Mapping plate_name → 'wx wy wz' string (DEG/My).
    igb14_folder   : str   Folder containing ``*_igb14.vel`` files.
    igb_nocomb_path: str   Sub-path segment (e.g. ``'igb14_no_comb'``).
    results_path   : str   Root results folder (e.g. ``'results'``).
    pycvframe_path : str   Absolute path to pycvframe.py.
    """
    dest_folders = {
        plate: os.path.join(results_path, igb_nocomb_path, plate)
        for plate in euler_poles
    }
    for folder in dest_folders.values():
        os.makedirs(folder, exist_ok=True)

    for vel_file in os.listdir(igb14_folder):
        if not vel_file.endswith('.vel'):
            continue
        base     = os.path.splitext(vel_file)[0]
        new_base = base.replace('_igb14', '')

        for plate, euler_params in euler_poles.items():
            out_name = f'{new_base}_{plate}.vel'
            cmd = [sys.executable, pycvframe_path, vel_file, out_name,
                   'ITRF14', euler_params]
            result = subprocess.run(
                cmd, cwd=igb14_folder,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
            out_path = os.path.join(igb14_folder, out_name)
            if not os.path.exists(out_path):
                print(f"WARNING: pycvframe failed for {vel_file} → {out_name}")
                print(f"  stdout: {result.stdout.decode().strip()}")
                print(f"  stderr: {result.stderr.decode().strip()}")
                continue
            shutil.move(out_path, os.path.join(dest_folders[plate], out_name))

    print("Euler-pole rotations completed.")


def run_velocity_combination(frames, combination_script_path,
                              igb_nocomb_path, results_path,
                              output_folder, mode='auto'):
    """Run the combine_vel.py script for each reference frame.

    Parameters
    ----------
    frames                  : list[str]  Plate names plus ``'igb14'``.
    combination_script_path : str        Path to combine_vel.py.
    igb_nocomb_path         : str        Sub-path segment (e.g. ``'igb14_no_comb'``).
    results_path            : str        Root results folder.
    output_folder           : str        Destination for combined velocity CSVs.
    mode                    : str        ``'auto'``, ``'parallel'``, or
                                         ``'sequential'``.
    """
    if not os.path.exists(combination_script_path):
        raise FileNotFoundError(
            f"combine_vel.py not found: {combination_script_path}"
        )

    def _combine_one(ref_frame):
        print(f"Combining GPS velocities in {ref_frame}-fixed reference frame ...")
        folder = os.path.join(results_path, igb_nocomb_path, ref_frame)
        subprocess.run(
            ['python', combination_script_path, folder, output_folder],
            check=True,
        )

    print("--------------------------------------------------------------")
    print("               Combining rotated velocity fields              ")
    print("--------------------------------------------------------------")
    print()

    if mode == 'sequential':
        print("Sequential execution.")
        for frame in frames:
            _combine_one(frame)

    elif mode == 'parallel':
        n_workers = min(os.cpu_count() or 1, len(frames))
        print(f"Parallel execution ({n_workers} workers).")
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as pool:
            list(pool.map(_combine_one, frames))

    elif mode == 'auto':
        mem_gb = get_available_memory_gb()
        print(f"Available memory: {mem_gb:.2f} GB")
        if mem_gb < 1.0:
            print("Low memory — running sequentially.")
            for frame in frames:
                _combine_one(frame)
        else:
            n_workers = min(os.cpu_count() or 1, len(frames))
            print(f"Parallel execution ({n_workers} workers).")
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as pool:
                list(pool.map(_combine_one, frames))

    else:
        raise ValueError(
            f"Invalid mode '{mode}'. Use 'auto', 'parallel', or 'sequential'."
        )
