"""ITRF14 alignment pipeline using pyvelrot.

Functions
---------
align_to_itrf14     Run pyvelrot per .vel file to produce *_igb14.vel outputs.
collect_igb14_files Copy the per-file IGB14 results into the shared igb14 folder.
"""

import os
import shutil
import subprocess
import sys

# Sentinel strings that bracket the aligned-velocity block in pyvelrot output
_TOP = "SYSTEM 1 Velocities transformed to SYSTEM 2"
_BOT = "SYSTEM 2 Velocities except those in SYSTEM 1"


def align_to_itrf14(input_vel_folder, rotation_folder, reference_vel_path,
                    lnk_file_path, pyvelrot_path):
    """Run pyvelrot for each .vel file to produce ITRF14-aligned outputs.

    For every ``*.vel`` file in *input_vel_folder*:
    - Creates a per-file subdirectory inside *rotation_folder*.
    - Copies the .vel file, the reference velocity file, and the link file.
    - If the file IS the reference, renames it directly to ``*_igb14.vel``.
    - Otherwise runs pyvelrot via subprocess, parses the tmp output, and
      writes the extracted lines to ``*_igb14.vel``.

    Parameters
    ----------
    input_vel_folder   : str  Folder containing ``*.vel`` files to align.
    rotation_folder    : str  Root output folder; per-file subdirs are created here.
    reference_vel_path : str  Path to the reference .vel file (ITRF14 anchor).
    lnk_file_path      : str  Path to the velrot.lnk link file.
    pyvelrot_path      : str  Absolute path to pyvelrot.py.
    """
    if not os.path.exists(reference_vel_path):
        raise FileNotFoundError(
            f"Reference velocity file not found: {reference_vel_path}"
        )

    if os.path.exists(rotation_folder):
        for entry in os.listdir(rotation_folder):
            ep = os.path.join(rotation_folder, entry)
            if os.path.isfile(ep):
                os.unlink(ep)
            elif os.path.isdir(ep):
                shutil.rmtree(ep)
    else:
        os.makedirs(rotation_folder)

    ref_name = os.path.basename(reference_vel_path)
    ref_base = os.path.splitext(ref_name)[0]
    lnk_name = os.path.basename(lnk_file_path)

    for vel_file in os.listdir(input_vel_folder):
        if not vel_file.endswith('.vel'):
            continue
        base = os.path.splitext(vel_file)[0]
        sub = os.path.join(rotation_folder, base)
        os.makedirs(sub, exist_ok=True)

        shutil.copy(os.path.join(input_vel_folder, vel_file), sub)
        shutil.copy(reference_vel_path, sub)
        shutil.copy(lnk_file_path, sub)

        igb14_vel = os.path.join(sub, f'{base}_igb14.vel')

        if base == ref_base:
            # Reference file is already in ITRF14; just rename it.
            os.rename(os.path.join(sub, ref_name), igb14_vel)
            continue

        tmp_name = f'{base}_igb14.tmp'
        cmd = [
            sys.executable, pyvelrot_path,
            vel_file, 'NONE', ref_name,
            'NONE', tmp_name, 'NONE', lnk_name, '0', 'TR',
        ]
        log_path = os.path.join(sub, f'{base}_align.log')
        with open(log_path, 'w') as log_fh:
            subprocess.run(cmd, cwd=sub, stdout=log_fh,
                           stderr=subprocess.STDOUT, check=True)

        tmp_path = os.path.join(sub, tmp_name)
        with open(tmp_path) as fh:
            lines = fh.readlines()

        start = next(i for i, l in enumerate(lines) if _TOP in l)
        end   = next(i for i, l in enumerate(lines) if _BOT in l)
        body  = [l.replace('*', '').rstrip() for l in lines[start + 3: end - 1]]

        with open(igb14_vel, 'w') as fh:
            fh.write('\n'.join(body))


def collect_igb14_files(rotation_folder, igb14_dest_folder):
    """Copy ``*_igb14.vel`` files from rotation sub-folders to a shared folder.

    Parameters
    ----------
    rotation_folder    : str  Root folder created by :func:`align_to_itrf14`.
    igb14_dest_folder  : str  Destination (e.g. ``results/igb14_no_comb/igb14``).
    """
    os.makedirs(igb14_dest_folder, exist_ok=True)
    for folder in os.listdir(rotation_folder):
        sub = os.path.join(rotation_folder, folder)
        if not os.path.isdir(sub):
            continue
        src = os.path.join(sub, f'{folder}_igb14.vel')
        if os.path.exists(src):
            shutil.copy(src, igb14_dest_folder)
