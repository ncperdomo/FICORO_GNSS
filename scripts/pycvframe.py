"""
PyCvframe - Python implementation of GAMIT cvframe
===================================================

Rotates a GNSS velocity field from one reference frame to another by
subtracting the velocity field predicted by a known Euler rotation pole.

This module reproduces the behavior of the GAMIT/GLOBK Fortran program
cvframe.f written by T. Herring et al.

Usage as a module:
    from pycvframe import PyCvframe
    cv = PyCvframe()
    cv.run("input.vel", "output.vel", rotation_pole=[wx, wy, wz])

Usage as a script:
    python pycvframe.py input.vel output.vel INFRAME "wx wy wz"
"""

import numpy as np
import sys
import os
from typing import Union, List, Tuple


# ---------------------------------------------------------------------------
# Constants (from GAMIT const_param.h)
# ---------------------------------------------------------------------------

PI = 3.1415926535897932
EARTH_RAD = 6378137.0          # WGS-84 semi-major axis (m)
EARTH_FLAT = 1.0 / 298.257222101  # WGS-84 flattening
EARTH_E2 = 2 * EARTH_FLAT - EARTH_FLAT ** 2  # first eccentricity squared


# ---------------------------------------------------------------------------
# Coordinate and rotation utilities
# ---------------------------------------------------------------------------

def geod_to_xyz(lat_deg: float, lon_deg: float, height: float = 0.0) -> np.ndarray:
    """
    Convert geodetic (lat, lon, height) to ECEF XYZ.

    Mirrors GAMIT's GEOD_to_XYZ subroutine (which receives co-latitude).
    Height is assumed zero (sites on the ellipsoid surface).
    """
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)
    N = EARTH_RAD / np.sqrt(1.0 - EARTH_E2 * sin_lat ** 2)
    X = (N + height) * cos_lat * np.cos(lon)
    Y = (N + height) * cos_lat * np.sin(lon)
    Z = (N * (1.0 - EARTH_E2) + height) * sin_lat
    return np.array([X, Y, Z])


def rotation_matrix_neu(lat_deg: float, lon_deg: float) -> np.ndarray:
    """
    Compute 3×3 rotation matrix R such that V_neu = R @ V_xyz.

    Row order: North (0), East (1), Up (2).
    This is the standard geodetic NEU convention used by GAMIT's rotate_geod.
    """
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    sin_lat, cos_lat = np.sin(lat), np.cos(lat)
    sin_lon, cos_lon = np.sin(lon), np.cos(lon)
    return np.array([
        [-sin_lat * cos_lon, -sin_lat * sin_lon,  cos_lat],  # North
        [-sin_lon,            cos_lon,             0.0    ],  # East
        [ cos_lat * cos_lon,  cos_lat * sin_lon,  sin_lat],  # Up
    ])


def cross_product(omega: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """Return omega × pos (velocity due to solid-body rotation)."""
    return np.cross(omega, pos)


def _batch_transform(lons: np.ndarray, lats: np.ndarray,
                     ve: np.ndarray, vn: np.ndarray,
                     rot_prime: np.ndarray):
    """
    Vectorized equivalent of apply_frame_rotation for N sites.

    All trig is computed once over the whole array; the cross product is
    inlined; only the North (R row 0) and East (R row 1) projections of
    the rotation matrix are evaluated — the Up component is never needed.
    """
    lat = np.radians(lats)
    lon = np.radians(lons)
    slat, clat = np.sin(lat), np.cos(lat)
    slon, clon = np.sin(lon), np.cos(lon)

    # geod_to_xyz — vectorized, (N,) each
    Nr = EARTH_RAD / np.sqrt(1.0 - EARTH_E2 * slat ** 2)
    px = Nr * clat * clon
    py = Nr * clat * slon
    pz = Nr * (1.0 - EARTH_E2) * slat

    # omega × pos — inline, avoids np.cross per-call overhead
    ox, oy, oz = rot_prime
    cvx = oy * pz - oz * py
    cvy = oz * px - ox * pz
    cvz = ox * py - oy * px

    # R[0] (North) = [-slat*clon, -slat*slon,  clat]
    # R[1] (East)  = [-slon,       clon,       0   ]
    neu_N = -slat * clon * cvx - slat * slon * cvy + clat * cvz
    neu_E = -slon * cvx        + clon * cvy

    return ve - neu_E * 1000.0, vn - neu_N * 1000.0


# ---------------------------------------------------------------------------
# Frame-to-frame rotation lookup — delegates to frame_registry
# ---------------------------------------------------------------------------

try:
    from frame_registry import frame_to_frame, list_frames, FRAME_REGISTRY, FRAME_NAMES
except ImportError:
    from .frame_registry import frame_to_frame, list_frames, FRAME_REGISTRY, FRAME_NAMES


# ---------------------------------------------------------------------------
# Core transformation
# ---------------------------------------------------------------------------

def apply_frame_rotation(lon: float, lat: float,
                         ve: float, vn: float,
                         rot_prime: np.ndarray) -> Tuple[float, float]:
    """
    Subtract the rotational velocity predicted by rot_prime (rad/yr) from
    the East and North velocity components at site (lon, lat).

    Parameters
    ----------
    lon, lat : float
        Site longitude and latitude in decimal degrees.
    ve, vn : float
        East and North velocity in mm/yr.
    rot_prime : np.ndarray, shape (3,)
        Euler rotation vector in rad/yr (XYZ components).

    Returns
    -------
    ve_new, vn_new : float
        Corrected East and North velocity in mm/yr.

    Notes
    -----
    Implements exactly the logic in cvframe.f:
        site_vel = cross_prod(rot_prime, site_pos)
        neu_vel  = rotate_geod(site_vel, XYZ -> NEU)
        ve -= neu_vel[1] * 1000   (East, index 2 in Fortran = index 1 in Python)
        vn -= neu_vel[0] * 1000   (North, index 1 in Fortran = index 0 in Python)
    The Up component is NOT modified.
    """
    site_pos = geod_to_xyz(lat, lon)
    site_vel = cross_product(rot_prime, site_pos)
    R = rotation_matrix_neu(lat, lon)
    neu_vel = R @ site_vel          # [N, E, U] in m/yr
    ve_new = ve - neu_vel[1] * 1000.0
    vn_new = vn - neu_vel[0] * 1000.0
    return ve_new, vn_new


# ---------------------------------------------------------------------------
# File I/O
# ---------------------------------------------------------------------------

def _parse_data_line(line: str):
    """
    Parse one data line from a .vel file.

    Returns a dict with keys: lon, lat, ve, vn, de, dn, se, sn, rho, vu, du, su, site
    or None if the line cannot be decoded.

    A data line starts with a space (first character).
    """
    tokens = line.split()
    if len(tokens) < 12:
        return None
    try:
        lon  = float(tokens[0])
        lat  = float(tokens[1])
        ve   = float(tokens[2])
        vn   = float(tokens[3])
        de   = float(tokens[4])
        dn   = float(tokens[5])
        se   = float(tokens[6])
        sn   = float(tokens[7])
        rho  = float(tokens[8])
        vu   = float(tokens[9])
        du   = float(tokens[10])
        su   = float(tokens[11])
        site = tokens[12] if len(tokens) > 12 else ""
    except (ValueError, IndexError):
        return None
    return dict(lon=lon, lat=lat, ve=ve, vn=vn, de=de, dn=dn,
                se=se, sn=sn, rho=rho, vu=vu, du=du, su=su, site=site)


def _format_data_line(lon, lat, ve, vn, de, dn, se, sn, rho, vu, du, su,
                      site: str) -> str:
    """
    Format a data line exactly as cvframe.f format 220:
        format(2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,3(1x,f7.2), 1x,a9)

    The outer 1x between lat and the 6-velocity group adds one extra space
    before ve; the 2x before the up-rate group adds one extra space before vu.
    """
    return (
        f" {lon:10.5f} {lat:10.5f}"
        f"  {ve:7.2f} {vn:7.2f} {de:7.2f} {dn:7.2f} {se:7.2f} {sn:7.2f}"
        f" {rho:6.3f}   {vu:7.2f} {du:7.2f} {su:7.2f} {site:<9s}"
    )


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class PyCvframe:
    """
    Python equivalent of GAMIT's cvframe program.

    Transforms a GNSS velocity field from one reference frame to another
    by subtracting the velocity predicted by a known Euler rotation pole.
    """

    def __init__(self):
        self.rot_prime = None       # rotation vector in rad/yr
        self.inframe = ""
        self.outframe = ""
        self.invel = ""
        self.outvel = ""

    def run(self,
            invel: str,
            outvel: str,
            inframe: str,
            outframe_or_pole: Union[str, List[float], np.ndarray],
            pole_units: str = "deg/Myr") -> None:
        """
        Transform velocity file invel → outvel.

        Parameters
        ----------
        invel : str
            Path to input .vel file.
        outvel : str
            Path to output .vel file.
        inframe : str
            Name of the input reference frame (e.g. "ITRF14", "NONE").
        outframe_or_pole : str or array-like
            Either the name of the target frame (looked up via frame_to_frame)
            or a 3-element Euler pole vector [wx, wy, wz].
        pole_units : str
            Units of the numeric pole, either "deg/Myr" (default, matching
            cvframe.f) or "rad/yr".
        """
        self.invel = invel
        self.outvel = outvel
        self.inframe = inframe.upper()

        # Resolve rotation pole
        if isinstance(outframe_or_pole, str):
            self.outframe = outframe_or_pole.upper()
            rot_raw = frame_to_frame(self.outframe, self.inframe)
            self.rot_prime = rot_raw  # already rad/yr from the lookup
        else:
            pole = np.asarray(outframe_or_pole, dtype=float)
            if pole_units == "deg/Myr":
                # Fortran: rot_prime = pole * pi/180 / 1e6
                self.rot_prime = pole * PI / 180.0 / 1.0e6
            else:
                self.rot_prime = pole
            self.outframe = "USER"

        # Fortran condition: abs(rot_prime(1))+abs(rot_prime(2))+abs(rot_prime(3)) < 1e-15
        if np.sum(np.abs(self.rot_prime)) < 1.0e-15:
            print(f"No rotation rate between {self.inframe} and {self.outframe}")
            return

        print(f"CVFRAME: Rotation Pole "
              f"{self.rot_prime[0]*180/PI*1e6:12.6f} "
              f"{self.rot_prime[1]*180/PI*1e6:12.6f} "
              f"{self.rot_prime[2]*180/PI*1e6:12.6f} deg/Myr")

        self._process(invel, outvel)

    def _process(self, invel: str, outvel: str) -> None:
        """Read input, transform, write output (mirrors cvframe.f main loop)."""
        with open(invel, "r") as fin, open(outvel, "w") as fout:
            # Write header lines (matching cvframe.f format statements 140/145/155)
            fout.write(
                f"* CVFRAME: Vel file {os.path.basename(invel)}"
                f" rotated from {self.inframe:<8s} to {self.outframe:<8s}\n"
            )
            fout.write(
                f"* Rotation Pole "
                f"{self.rot_prime[0]*180/PI*1e6:12.6f}"
                f"{self.rot_prime[1]*180/PI*1e6:12.6f}"
                f"{self.rot_prime[2]*180/PI*1e6:12.6f} deg/Myr\n"
            )
            fout.write(
                "*  Long.       Lat.        E & N Rate     "
                " E & N Adj.    E & N +-  RHO       "
                "H Rate  H adj.   +- SITE\n"
            )
            fout.write(
                "*  (deg)      (deg)         (mm/yr)      "
                "(mm/yr)      (mm/yr)               (mm/yr)\n"
            )

            # Process lines
            n_transformed = 0
            n_error = 0
            for line in fin:
                line = line.rstrip("\n")
                if line and line[0] == " ":
                    rec = _parse_data_line(line)
                    if rec is None:
                        n_error += 1
                        fout.write(line + "\n")
                        continue
                    ve_new, vn_new = apply_frame_rotation(
                        rec["lon"], rec["lat"],
                        rec["ve"], rec["vn"],
                        self.rot_prime
                    )
                    fout.write(
                        _format_data_line(
                            rec["lon"], rec["lat"],
                            ve_new, vn_new,
                            rec["de"], rec["dn"],
                            rec["se"], rec["sn"], rec["rho"],
                            rec["vu"], rec["du"], rec["su"],
                            rec["site"]
                        ) + "\n"
                    )
                    n_transformed += 1
                else:
                    # Pass through comment / header lines unchanged
                    # trimlen equivalent: write non-empty lines
                    fout.write((line or "") + "\n")

        print(f"  Transformed {n_transformed} sites -> {outvel}")
        if n_error:
            print(f"  {n_error} lines had decode errors and were passed through")

    def transform_array(self,
                        lons: np.ndarray,
                        lats: np.ndarray,
                        ve: np.ndarray,
                        vn: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Apply the frame rotation to arrays of velocities.

        Parameters
        ----------
        lons, lats : array-like
            Site longitudes and latitudes in decimal degrees.
        ve, vn : array-like
            East and North velocities in mm/yr.

        Returns
        -------
        ve_new, vn_new : np.ndarray
            Corrected velocities in mm/yr.
        """
        lons = np.asarray(lons, dtype=float)
        lats = np.asarray(lats, dtype=float)
        ve   = np.asarray(ve,   dtype=float)
        vn   = np.asarray(vn,   dtype=float)
        if self.rot_prime is None:
            raise RuntimeError(
                "No rotation pole set. Call run() or set rot_prime before transform_array()."
            )
        return _batch_transform(lons, lats, ve, vn, self.rot_prime)


# ---------------------------------------------------------------------------
# Command-line interface (mirrors cvframe runstring)
# ---------------------------------------------------------------------------

def main():
    """
    Command-line interface:
        python pycvframe.py <invel> <outvel> <inframe> <outframe_or_pole>

    where <outframe_or_pole> is either a frame name or a quoted string of
    three space-separated values in deg/Myr, e.g. "0.328 -0.035 0.407".
    """
    if len(sys.argv) < 5:
        print("Usage: pycvframe.py <invel> <outvel> <inframe> <outframe_or_pole>")
        print("  outframe_or_pole: frame name OR 'wx wy wz' in deg/Myr")
        sys.exit(1)

    invel  = sys.argv[1]
    outvel = sys.argv[2]
    inframe = sys.argv[3]
    pole_arg = sys.argv[4]

    # Decide if numeric pole or frame name (same logic as cvframe.f)
    parts = pole_arg.replace(":", " ").split()
    try:
        pole = [float(p) for p in parts]
        outframe_or_pole = np.array(pole)
        pole_units = "deg/Myr"
    except ValueError:
        outframe_or_pole = pole_arg
        pole_units = "deg/Myr"

    cv = PyCvframe()
    cv.run(invel, outvel, inframe, outframe_or_pole, pole_units=pole_units)


if __name__ == "__main__":
    main()
