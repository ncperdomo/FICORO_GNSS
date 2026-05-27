"""
PyVelrot - Python implementation of GAMIT velrot
=================================================

Compares two GNSS velocity fields and estimates the Helmert transformation
(translations + rotations, optionally + scale) between them using weighted
least squares.  The transformed System 1 velocities are written to an
output .vel file together with the System 2 sites.

This module reproduces the behavior of the GAMIT/GLOBK Fortran program
velrot.f written by T. Herring et al.

Usage as a module:
    from pyvelrot import PyVelrot
    vr = PyVelrot()
    vr.run("sys1.vel", "NONE", "sys2.vel", "NONE",
           "output.vel", "NONE", "links.lnk",
           height_weight=0, param_opt="TR")

Usage as a script:
    python pyvelrot.py sys1.vel NONE sys2.vel NONE out.vel NONE links.lnk 0 TR
"""

import numpy as np
import sys
import os
import datetime
from typing import List, Tuple, Optional, Dict


# ---------------------------------------------------------------------------
# Constants (from GAMIT const_param.h)
# ---------------------------------------------------------------------------

PI          = 3.1415926535897932
EARTH_RAD   = 6378137.0                     # WGS-84 semi-major axis (m)
EARTH_FLAT  = 1.0 / 298.257222101           # WGS-84 flattening
EARTH_E2    = 2 * EARTH_FLAT - EARTH_FLAT**2
RAD_TO_MAS  = 648000.0e3 / PI              # radians → milliarcseconds
VELROT_VER  = "1.01"

# Minimum height weight to avoid numerical problems (from velrot.f line 183)
MIN_HEIGHT_WEIGHT = 1.0e-6

# Tolerance for deciding to include vertical component (increment_norm)
HEIGHT_WEIGHT_TOL = 1.0e-5

# Large value used to constrain unused parameters (clear_norm)
HUGE_DIAGONAL = 1.0e14


# ---------------------------------------------------------------------------
# Coordinate utilities
# ---------------------------------------------------------------------------

def geod_to_xyz(lat_deg: float, lon_deg: float, height: float = 0.0) -> np.ndarray:
    """
    Geodetic (lat °, lon °, height m) → ECEF XYZ (m).
    Matches GAMIT GEOD_to_XYZ, which receives co-latitude internally.
    """
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)
    N = EARTH_RAD / np.sqrt(1.0 - EARTH_E2 * sin_lat**2)
    X = (N + height) * cos_lat * np.cos(lon)
    Y = (N + height) * cos_lat * np.sin(lon)
    Z = (N * (1.0 - EARTH_E2) + height) * sin_lat
    return np.array([X, Y, Z])


def xyz_to_geod(xyz: np.ndarray) -> Tuple[float, float, float]:
    """ECEF XYZ → (lat °, lon °, height m) via iterative Bowring."""
    X, Y, Z = xyz
    lon = np.degrees(np.arctan2(Y, X))
    p = np.sqrt(X**2 + Y**2)
    lat = np.arctan2(Z, p * (1.0 - EARTH_E2))
    for _ in range(5):
        sin_lat = np.sin(lat)
        N = EARTH_RAD / np.sqrt(1.0 - EARTH_E2 * sin_lat**2)
        height = p / np.cos(lat) - N
        lat = np.arctan2(Z, p * (1.0 - EARTH_E2 * N / (N + height)))
    return np.degrees(lat), lon, height


def rotation_matrix_neu(lat_deg: float, lon_deg: float) -> np.ndarray:
    """
    3×3 rotation matrix R such that V_NEU = R @ V_XYZ.
    Row 0 = North, Row 1 = East, Row 2 = Up.
    Matches GAMIT's rotate_geod (which uses co-latitude internally).
    """
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    sl, cl = np.sin(lat), np.cos(lat)
    sn, cn = np.sin(lon), np.cos(lon)
    return np.array([
        [-sl * cn, -sl * sn,  cl],   # North
        [-sn,       cn,       0.0],  # East
        [ cl * cn,  cl * sn,  sl],   # Up
    ])


def rotate_geod(vec_in: np.ndarray, from_frame: str, to_frame: str,
                pos_xyz: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Rotate between XYZ and NEU, also returning the local geodetic coordinates
    and the rotation matrix.

    Returns (vec_out, loc_coord, rot_matrix) where:
        loc_coord = (co_lat_rad, lon_rad, height_m)  — matches GAMIT convention
        rot_matrix = 3×3 matrix (R: XYZ→NEU)
    """
    lat_deg, lon_deg, height = xyz_to_geod(pos_xyz)
    R = rotation_matrix_neu(lat_deg, lon_deg)
    if from_frame.upper() == "XYZ" and to_frame.upper() == "NEU":
        vec_out = R @ vec_in
    elif from_frame.upper() == "NEU" and to_frame.upper() == "XYZ":
        vec_out = R.T @ vec_in
    else:
        raise ValueError(f"Invalid frame pair: {from_frame} -> {to_frame}")
    co_lat = np.radians(90.0 - lat_deg)
    lon_rad = np.radians(lon_deg)
    loc_coord = np.array([co_lat, lon_rad, height])
    return vec_out, loc_coord, R


def cross_product(omega: np.ndarray, pos: np.ndarray) -> Tuple[np.ndarray, float]:
    """omega × pos; also returns |sin(angle)| for diagnostics."""
    vel = np.cross(omega, pos)
    mag_o = np.linalg.norm(omega)
    mag_p = np.linalg.norm(pos)
    sang = np.linalg.norm(vel) / (mag_o * mag_p) if (mag_o > 0 and mag_p > 0) else 0.0
    return vel, sang


# ---------------------------------------------------------------------------
# Frame lookup — delegates to frame_registry
# ---------------------------------------------------------------------------

try:
    from frame_registry import frame_to_frame, list_frames, FRAME_REGISTRY, FRAME_NAMES
except ImportError:
    from .frame_registry import frame_to_frame, list_frames, FRAME_REGISTRY, FRAME_NAMES


# ---------------------------------------------------------------------------
# Velocity file I/O
# ---------------------------------------------------------------------------

def read_vel_file(filename: str) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Read a GAMIT .vel file.

    Returns
    -------
    coords : np.ndarray, shape (N, 3, 2)
        coords[i, :, 0] = XYZ position  (m)
        coords[i, :, 1] = XYZ velocity  (m/yr)
    covs : np.ndarray, shape (N, 3, 3)
        NEU velocity covariance matrix in (m/yr)²:
            cov[i, 0, 0] = var_N,  cov[i, 1, 1] = var_E,  cov[i, 2, 2] = var_U
            cov[i, 0, 1] = cov[i, 1, 0] = cov_NE
    names : list of str (length N)

    Notes
    -----
    Mirrors velrot.f subroutine read_sys.
    Data lines start with a space; any other first character is a comment.
    Column layout: lon lat ve vn de dn se sn rho vu du su name
    Velocities are stored internally in m/yr, covariances in (m/yr)².
    """
    coords_list = []
    covs_list = []
    names_list = []

    with open(filename, "r") as fh:
        for line in fh:
            if not line or line[0] != " ":
                continue
            line = line.strip()
            if not line:
                continue
            tokens = line.split()
            if len(tokens) < 12:
                continue
            try:
                vals = [float(t) for t in tokens[:12]]
            except ValueError:
                continue

            lon, lat = vals[0], vals[1]
            ve, vn, de, dn = vals[2], vals[3], vals[4], vals[5]
            se, sn, rho    = vals[6], vals[7], vals[8]
            vu, du, su     = vals[9], vals[10], vals[11]
            site = tokens[12] if len(tokens) > 12 else ""

            # Convert geodetic position to XYZ (height = 0)
            pos_xyz = geod_to_xyz(lat, lon)

            # Velocity in NEU (m/yr): NEU order = [N, E, U]
            vel_neu = np.array([vn / 1000.0, ve / 1000.0, vu / 1000.0])

            # Rotate to XYZ velocity
            vel_xyz, _, R = rotate_geod(vel_neu, "NEU", "XYZ", pos_xyz)

            # NEU covariance in (m/yr)²  (index 0=N, 1=E, 2=U)
            cov = np.zeros((3, 3))
            cov[0, 0] = sn**2 * 1.0e-6        # var_N
            cov[1, 1] = se**2 * 1.0e-6        # var_E
            cov[0, 1] = rho * sn * se * 1.0e-6
            cov[1, 0] = cov[0, 1]
            cov[2, 2] = su**2 * 1.0e-6        # var_U

            coord = np.zeros((3, 2))
            coord[:, 0] = pos_xyz
            coord[:, 1] = vel_xyz

            coords_list.append(coord)
            covs_list.append(cov)
            names_list.append(site)

    if not coords_list:
        return np.zeros((0, 3, 2)), np.zeros((0, 3, 3)), []

    return (np.array(coords_list),
            np.array(covs_list),
            names_list)


# ---------------------------------------------------------------------------
# get_parts: design matrix
# ---------------------------------------------------------------------------

def get_parts(coord_xyz: np.ndarray, rot_matrix: np.ndarray,
              num_parn: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build the 3×7 design matrix for the Helmert transformation.

    Parameters
    ----------
    coord_xyz : (3,) position in XYZ (m)
    rot_matrix : (3, 3) rotation XYZ→NEU
    num_parn : number of parameters (3, 6, or 7)

    Returns
    -------
    A_neu : (3, 7) in NEU frame
    A_xyz : (3, 7) in XYZ frame

    Columns: [Tx, Ty, Tz, Rx, Ry, Rz, scale]

    Notes
    -----
    Mirrors velrot.f subroutine get_parts.
    Rotation partials use the linearised convention:
        dV/dRx = [0, -Z,  Y]^T
        dV/dRy = [Z,  0, -X]^T
        dV/dRz = [-Y, X,  0]^T
    (i.e. V += omega × X  with the sign convention of the Fortran code).
    """
    X, Y, Z = coord_xyz
    A_xyz = np.zeros((3, 7))

    # Translations
    A_xyz[0, 0] = 1.0
    A_xyz[1, 1] = 1.0
    A_xyz[2, 2] = 1.0

    # Rotations (omega × X)
    A_xyz[1, 3] = -Z
    A_xyz[2, 3] =  Y
    A_xyz[0, 4] =  Z
    A_xyz[2, 4] = -X
    A_xyz[0, 5] = -Y
    A_xyz[1, 5] =  X

    # Scale
    if num_parn == 7:
        A_xyz[0, 6] = X
        A_xyz[1, 6] = Y
        A_xyz[2, 6] = Z

    # Rotate to NEU
    A_neu = rot_matrix @ A_xyz
    A_neu[:, num_parn:] = 0.0   # zero out columns beyond num_parn
    return A_neu, A_xyz


# ---------------------------------------------------------------------------
# var_comp: variance propagation  A C A^T
# ---------------------------------------------------------------------------

def var_comp(A_neu: np.ndarray, C: np.ndarray) -> np.ndarray:
    """
    Propagate covariance: result = A_neu @ C @ A_neu.T
    Mirrors GAMIT's var_comp (kinv_lib) called with iopt=1.
    """
    return A_neu @ C @ A_neu.T


# ---------------------------------------------------------------------------
# PyVelrot class
# ---------------------------------------------------------------------------

class PyVelrot:
    """
    Python equivalent of GAMIT velrot.

    State variables mirror the Fortran COMMON block /corcom_com/.
    """

    def __init__(self):
        # System data (set by read_vel_file)
        self.sys1_coord: np.ndarray = None  # (N1, 3, 2)
        self.sys1_cov:   np.ndarray = None  # (N1, 3, 3)
        self.sys1_names: List[str]  = []
        self.sys2_coord: np.ndarray = None  # (N2, 3, 2)
        self.sys2_cov:   np.ndarray = None  # (N2, 3, 3)
        self.sys2_names: List[str]  = []

        # Transformation settings
        self.param_opt: str   = "TR"
        self.num_parn:  int   = 6
        self.height_weight: float = 1.0

        # Frame strings
        self.sys1_frame: str = "NONE"
        self.sys2_frame: str = "NONE"
        self.out_frame:  str = "NONE"

        # Fundamental links [(i_sys1, i_sys2), ...]
        self.fund_links: List[Tuple[int, int]] = []
        self.eq_dist: float = 0.0
        self.cp_dist: float = 0.0
        self.av_dist: float = 0.0

        # Estimation results
        self.norm_eq:       np.ndarray = np.zeros((7, 7))
        self.bvec:          np.ndarray = np.zeros(7)
        self.trans_parm:    np.ndarray = np.zeros(7)   # internal units
        self.trans_parm_out: np.ndarray = np.zeros(7)  # display units
        self.trans_sigma:   np.ndarray = np.zeros(7)   # display units

        self.sum_prefit:  float = 0.0
        self.sum_postfit: float = 0.0
        self.sum_weight:  float = 0.0
        self.chi_fit:     float = 1.0
        self.rms_fit:     float = 0.0
        self.num_data:    int   = 0

        # Output
        self.out_file: str  = ""
        self.out_unit       = None   # file handle or None (→ stdout)

    # -----------------------------------------------------------------------
    # Main entry point
    # -----------------------------------------------------------------------

    def run(self,
            sys1_file: str,
            sys1_frame: str,
            sys2_file: str,
            sys2_frame: str,
            out_file: str,
            out_frame: str,
            fund_file: str = "",
            height_weight: float = 1.0,
            param_opt: str = "TR") -> None:
        """
        Full velrot run.

        Parameters
        ----------
        sys1_file : str
            Path to System 1 .vel file.
        sys1_frame : str
            Reference frame of System 1 (e.g. "NONE", "ITRF14").
        sys2_file : str
            Path to System 2 .vel file (the reference).
        sys2_frame : str
            Reference frame of System 2.
        out_file : str
            Path to the output .vel file.
        out_frame : str
            Frame for output velocities.
        fund_file : str
            Path to the fundamental sites file (or "" to use all common sites).
        height_weight : float
            Weight for the vertical component. 0 = no vertical, 1 = equal.
            Values < 1e-6 are clamped to 1e-6 (Fortran behaviour).
        param_opt : str
            Transformation type: "T" (3-param), "TR" (6-param, default),
            "TRS" (7-param), "L" (2-param local).
        """
        # --- Parse param_opt --------------------------------------------------
        self.param_opt  = param_opt.upper()
        self.num_parn   = self._parse_param_opt(self.param_opt)
        self.sys1_frame = sys1_frame.strip().upper()
        self.sys2_frame = sys2_frame.strip().upper()
        self.out_frame  = out_frame.strip().upper() if out_frame.strip() else self.sys1_frame
        self.out_file   = out_file

        # Clamp height_weight (mirrors velrot.f line 183)
        if height_weight < MIN_HEIGHT_WEIGHT:
            height_weight = MIN_HEIGHT_WEIGHT
        self.height_weight = height_weight

        print(f"\n VELROT: Velocity field comparison and combination Version {VELROT_VER}\n")

        # --- Read sys1 --------------------------------------------------------
        self.sys1_coord, self.sys1_cov, self.sys1_names = read_vel_file(sys1_file)
        print(f" There are {len(self.sys1_names):5d} sites in sys file {sys1_file}")

        # --- Frame update sys1 → out_frame ------------------------------------
        frame_epoch = 2000.0
        self._frame_update(1, self.sys1_frame, self.out_frame, frame_epoch)

        # --- Read sys2 --------------------------------------------------------
        self.sys2_coord, self.sys2_cov, self.sys2_names = read_vel_file(sys2_file)
        print(f" There are {len(self.sys2_names):5d} sites in sys file {sys2_file}")

        # --- Frame update sys2 → out_frame ------------------------------------
        self._frame_update(2, self.sys2_frame, self.out_frame, frame_epoch)

        # --- Read fundamental sites -------------------------------------------
        self._read_fund_sites(fund_file)

        # --- Estimate transformation ------------------------------------------
        self._transframe()

        # --- Output summary ---------------------------------------------------
        fh = self._open_out(out_file)
        self._output_sum(fh, sys1_file, sys2_file, fund_file, frame_epoch)

        # --- Apply transformation to sys1 ------------------------------------
        self._update_tran()

        # --- Write transformed velocity field --------------------------------
        self._output_frame(fh)

        if fh is not sys.stdout:
            fh.close()

    # -----------------------------------------------------------------------
    # Parameter parsing
    # -----------------------------------------------------------------------

    @staticmethod
    def _parse_param_opt(opt: str) -> int:
        n = 0
        if "T" in opt: n = 3
        if "R" in opt: n = 6
        if "S" in opt: n = 7
        if "L" in opt: n = 2
        if n == 0:
            n = 6
        return n

    # -----------------------------------------------------------------------
    # Frame update
    # -----------------------------------------------------------------------

    def _frame_update(self, which: int, sys_frame: str, out_frm: str,
                      frame_epoch: float) -> None:
        """
        Rotate system velocities into out_frame.
        Mirrors velrot.f subroutine frame_update.
        """
        rot_vec = frame_to_frame(sys_frame, out_frm)

        print(f" Rotating from {sys_frame:<10s} to {out_frm:<10s} using rotation"
              f" vector {rot_vec[0]*180e6/PI:12.6f} {rot_vec[1]*180e6/PI:12.6f}"
              f" {rot_vec[2]*180e6/PI:12.6f} degs/Myrs")

        if np.sum(np.abs(rot_vec)) < 1.0e-10:
            return

        if which == 1:
            coords = self.sys1_coord
        else:
            coords = self.sys2_coord

        for i in range(len(coords)):
            xyz_vel, _ = cross_product(rot_vec, coords[i, :, 0])
            coords[i, :, 1] += xyz_vel

    # -----------------------------------------------------------------------
    # Fundamental sites
    # -----------------------------------------------------------------------

    def _read_fund_sites(self, fund_file: str) -> None:
        """
        Build fund_links list and set eq_dist/cp_dist/av_dist.
        Mirrors velrot.f subroutine read_fund_sites.
        """
        n1 = len(self.sys1_names)
        n2 = len(self.sys2_names)

        if fund_file and os.path.isfile(fund_file):
            links = []
            with open(fund_file, "r") as fh:
                for line in fh:
                    if not line or line[0] != " ":
                        continue
                    tokens = line.strip().split()
                    if not tokens:
                        continue
                    kw = tokens[0].upper()

                    if kw == "EQ_DIST":
                        if len(tokens) >= 2:
                            self.eq_dist = float(tokens[1])
                            self.cp_dist = self.eq_dist
                        # Build links for all pairs within eq_dist
                        for j in range(n1):
                            for k in range(n2):
                                dl = np.linalg.norm(
                                    self.sys1_coord[j, :, 0] - self.sys2_coord[k, :, 0]
                                )
                                if dl < self.eq_dist:
                                    links.append((j, k))

                    elif kw == "NAMES":
                        # Match all sites by name
                        links = []
                        name2_idx = {n: i for i, n in enumerate(self.sys2_names)}
                        for j, n1name in enumerate(self.sys1_names):
                            if n1name in name2_idx:
                                links.append((j, name2_idx[n1name]))

                    elif kw == "CP_DIST":
                        if len(tokens) >= 2:
                            self.cp_dist = float(tokens[1])

                    elif kw == "AV_DIST":
                        if len(tokens) >= 2:
                            self.av_dist = float(tokens[1])

                    else:
                        # Named site pair (or single name)
                        name1 = tokens[0]
                        name2 = tokens[1] if len(tokens) >= 2 else name1

                        # Handle +/- prefixes
                        remove1 = name1.startswith("-")
                        if remove1:
                            name1 = name1[1:]
                        elif name1.startswith("+"):
                            name1 = name1[1:]

                        remove2 = name2.startswith("-")
                        if remove2:
                            name2 = name2[1:]
                        elif name2.startswith("+"):
                            name2 = name2[1:]

                        # Apply case-folding to uppercase
                        name1 = name1.upper()
                        name2 = name2.upper()

                        ns1 = next(
                            (i for i, n in enumerate(self.sys1_names)
                             if n.upper() == name1), -1)
                        ns2 = next(
                            (i for i, n in enumerate(self.sys2_names)
                             if n.upper() == name2), -1)

                        if remove1 and ns1 >= 0:
                            links = [(a, b) for a, b in links if a != ns1]
                        if remove2 and ns2 >= 0:
                            links = [(a, b) for a, b in links if b != ns2]
                        if not remove1 and not remove2 and ns1 >= 0 and ns2 >= 0:
                            links.append((ns1, ns2))

            # Remove duplicates (preserve order)
            seen = set()
            unique = []
            for pair in links:
                key = (min(pair), max(pair))
                if key not in seen:
                    seen.add(key)
                    unique.append(pair)
            self.fund_links = unique

        else:
            # No file: match all by name
            name2_idx = {n: i for i, n in enumerate(self.sys2_names)}
            self.fund_links = [
                (j, name2_idx[n])
                for j, n in enumerate(self.sys1_names)
                if n in name2_idx
            ]

        print(f" There are {len(self.fund_links):5d} matching sites in fundamental"
              f" file {fund_file}")

    # -----------------------------------------------------------------------
    # Transformation estimation
    # -----------------------------------------------------------------------

    def _clear_norm(self) -> None:
        """Initialise normal equations with constraints on unused parameters."""
        self.norm_eq = np.zeros((7, 7))
        self.bvec    = np.zeros(7)
        self.sum_prefit  = 0.0
        self.sum_weight  = 0.0
        self.num_data    = 0

        # Constrain un-estimated parameters with large diagonal
        if "T" not in self.param_opt and "L" not in self.param_opt:
            for i in range(3):
                self.norm_eq[i, i] = HUGE_DIAGONAL
        if "S" not in self.param_opt:
            self.norm_eq[6, 6] = HUGE_DIAGONAL

    def _increment_norm(self, dn: np.ndarray, A_neu: np.ndarray,
                        weights: np.ndarray) -> None:
        """
        Accumulate normal equations for one site pair.
        Mirrors velrot.f subroutine increment_norm.

        Uses 2 components if w_U / w_E < HEIGHT_WEIGHT_TOL, else 3.
        """
        # Decide how many components to use
        if abs(weights[2] / weights[1]) < HEIGHT_WEIGHT_TOL:
            num_use = 2
        else:
            num_use = 3

        for i in range(num_use):
            for j in range(7):
                self.bvec[j]         += A_neu[i, j] * dn[i] * weights[i]
                for k in range(7):
                    self.norm_eq[j, k] += A_neu[i, j] * weights[i] * A_neu[i, k]
            self.sum_prefit += dn[i]**2 * weights[i]
            self.sum_weight += weights[i]
            self.num_data   += 1

    def _transframe(self) -> None:
        """
        Estimate transformation parameters by weighted least squares.
        Mirrors velrot.f subroutine transframe.
        After this call, self.norm_eq contains N^{-1} (the inverted matrix).
        """
        self._clear_norm()

        for ns1, ns2 in self.fund_links:
            # Velocity difference in XYZ
            dx = self.sys2_coord[ns2, :, 1] - self.sys1_coord[ns1, :, 1]

            # Rotate to NEU at sys1 position
            dn, loc_coord, R = rotate_geod(dx, "XYZ", "NEU",
                                           self.sys1_coord[ns1, :, 0])

            # Design matrix
            A_neu, _ = get_parts(self.sys1_coord[ns1, :, 0], R, self.num_parn)

            # Weights: inverse of combined NEU variance
            w = np.array([
                1.0 / (self.sys1_cov[ns1, 0, 0] + self.sys2_cov[ns2, 0, 0]),
                1.0 / (self.sys1_cov[ns1, 1, 1] + self.sys2_cov[ns2, 1, 1]),
                self.height_weight /
                    (self.sys1_cov[ns1, 2, 2] + self.sys2_cov[ns2, 2, 2]),
            ])

            self._increment_norm(dn, A_neu, w)

        # --- Solve normal equations ------------------------------------------
        b_save = self.bvec.copy()
        N_save = self.norm_eq.copy()

        # Solve and invert in one step using numpy
        # (mirrors invert_vis which both solves N*x=b and replaces N with N^-1)
        try:
            self.trans_parm = np.linalg.solve(N_save, b_save)
            self.norm_eq    = np.linalg.inv(N_save)
        except np.linalg.LinAlgError:
            self.trans_parm = np.zeros(7)
            self.norm_eq    = np.zeros((7, 7))

        # --- Post-fit statistics ---------------------------------------------
        # dprefit = b^T * N^{-1} * b  (= b^T * x)
        dprefit = b_save @ self.norm_eq @ b_save
        self.sum_postfit = self.sum_prefit - dprefit

        dof = self.num_data - self.num_parn
        if dof > 0:
            self.chi_fit = np.sqrt(self.sum_postfit / dof)
            self.rms_fit = np.sqrt(self.num_data / self.sum_weight) * self.chi_fit
            # Clamp chi to >= 0.1 (velrot.f line ~816)
            if self.chi_fit < 0.10:
                print(f"# Chi of fit less than 0.10 ({self.chi_fit:.4f}). Resetting to 1.0")
                self.chi_fit = 1.0
        else:
            self.chi_fit = 1.0
            self.rms_fit = 0.0
            self.trans_parm[:] = 0.0
            self.norm_eq[:] = 0.0

        # --- Convert to output units -----------------------------------------
        self.trans_parm_out    = np.zeros(7)
        self.trans_parm_out[0] = self.trans_parm[0] * 1000.0
        self.trans_parm_out[1] = self.trans_parm[1] * 1000.0
        self.trans_parm_out[2] = self.trans_parm[2] * 1000.0
        self.trans_parm_out[3] = self.trans_parm[3] * RAD_TO_MAS
        self.trans_parm_out[4] = self.trans_parm[4] * RAD_TO_MAS
        self.trans_parm_out[5] = self.trans_parm[5] * RAD_TO_MAS
        self.trans_parm_out[6] = self.trans_parm[6] * 1.0e9

        # Sigmas: sqrt of diagonal of N^{-1} scaled by chi_fit
        self.trans_sigma    = np.zeros(7)
        self.trans_sigma[0] = np.sqrt(max(0, self.norm_eq[0, 0])) * 1000.0      * self.chi_fit
        self.trans_sigma[1] = np.sqrt(max(0, self.norm_eq[1, 1])) * 1000.0      * self.chi_fit
        self.trans_sigma[2] = np.sqrt(max(0, self.norm_eq[2, 2])) * 1000.0      * self.chi_fit
        self.trans_sigma[3] = np.sqrt(max(0, self.norm_eq[3, 3])) * RAD_TO_MAS  * self.chi_fit
        self.trans_sigma[4] = np.sqrt(max(0, self.norm_eq[4, 4])) * RAD_TO_MAS  * self.chi_fit
        self.trans_sigma[5] = np.sqrt(max(0, self.norm_eq[5, 5])) * RAD_TO_MAS  * self.chi_fit
        self.trans_sigma[6] = np.sqrt(max(0, self.norm_eq[6, 6])) * 1.0e9       * self.chi_fit

        # Convert rms to mm/yr
        self.rms_fit *= 1000.0

    # -----------------------------------------------------------------------
    # Output summary
    # -----------------------------------------------------------------------

    _PARAM_LABS = ["X-Offset", "Y-Offset", "Z-Offset",
                   "X-Rot   ", "Y-Rot   ", "Z-Rot   ", "Scale   "]
    _UNIT_LABS  = ["(mm/yr)", "(mm/yr)", "(mm/yr)",
                   "(mas/yr)", "(mas/yr)", "(mas/yr)", "(ppb/yr)"]
    _COMP_LABS  = ["North", "East", "Up", "Horz"]

    def _open_out(self, out_file: str):
        if out_file.strip() in ("", "6"):
            return sys.stdout
        fh = open(out_file, "w")
        return fh

    def _output_sum(self, fh, sys1_file, sys2_file, fund_file, frame_epoch):
        """
        Write the velrot summary header and per-site residuals.
        Mirrors velrot.f subroutine output_sum.
        """
        now = datetime.datetime.now()
        fh.write(
            f"* VELROT Run on {now.year:4d}/{now.month:2d}/{now.day:2d}"
            f" {now.hour:2d}:{now.minute:2d} Version {VELROT_VER}\n"
        )
        fh.write(f"* SYSTEM 1 File    : {self.sys1_frame:<10s} {sys1_file}\n")
        fh.write(f"* SYSTEM 2 File    : {self.sys2_frame:<10s} {sys2_file}\n")
        fh.write(f"* FUNDAMENTAL File : {' ':11s}{fund_file}\n")
        fh.write(f"* OUTPUT FRAME     : {self.out_frame:<10s}"
                 f"  PARAM_OPT {self.param_opt}\n")
        fh.write(f"* EQ_DIST          :      {self.eq_dist:12.1f} m,"
                 f" CP_DIST : {self.cp_dist:12.1f} m*\n")
        fh.write(f"* AV_DIST          :      {self.av_dist:12.1f} m"
                 f"* HEIGHT WEIGHT    : {self.height_weight:10.6f}\n")
        fh.write("* \n")

        if self.num_data > 0:
            fh.write(
                f"* RMS fit for {self.num_data:8d} components"
                f" from {len(self.fund_links):4d} stations was"
                f" {self.rms_fit:10.2f} mm/yr, NRMS {self.chi_fit:8.2f}\n"
            )
            fh.write("* Estimates of Transformation parameters are: \n")
            for i in range(self.num_parn):
                fh.write(
                    f"*{i+1:5d} {self._PARAM_LABS[i]}"
                    f"   {self.trans_parm_out[i]:10.4f}"
                    f" {self.trans_sigma[i]:10.4f}"
                    f" {self._UNIT_LABS[i]}\n"
                )

        # Per-site residuals
        fh.write(
            "*  Differences at the fundamental sites\n"
            "*   #  Name 1  Name Ref    dN (mm)    dE (mm)"
            "    dU (mm)    sN (mm)    sE (mm)    sU (mm)"
            "   sTN (mm)   sTE (mm)   sTU (mm)\n"
        )

        # Statistics accumulators (mirrors velrot.f output_sum)
        summ = np.zeros(4)
        sumv = np.zeros(4)
        sumw = np.zeros(4)

        for idx, (ns1, ns2) in enumerate(self.fund_links):
            dx = self.sys2_coord[ns2, :, 1] - self.sys1_coord[ns1, :, 1]
            dn, loc, R = rotate_geod(dx, "XYZ", "NEU",
                                     self.sys1_coord[ns1, :, 0])

            A_neu, A_xyz = get_parts(self.sys1_coord[ns1, :, 0], R, self.num_parn)

            # Apply transformation: subtract predicted correction in XYZ
            for k in range(self.num_parn):
                dx -= A_xyz[:, k] * self.trans_parm[k]

            dn, _, _ = rotate_geod(dx, "XYZ", "NEU", self.sys1_coord[ns1, :, 0])

            # Propagated transformation uncertainty
            cov_neu = var_comp(A_neu, self.norm_eq)  # (3,3) in (m/yr)²

            # Combined sigma: sqrt(sys1_var + sys2_var)
            s_comb = np.array([
                np.sqrt(self.sys1_cov[ns1, j, j] + self.sys2_cov[ns2, j, j]) * 1000.0
                for j in range(3)
            ])
            s_tran = np.array([np.sqrt(max(0, cov_neu[j, j])) * 1000.0
                               for j in range(3)])

            fh.write(
                f"A{idx+1:5d} {self.sys1_names[ns1]:<8s} {self.sys2_names[ns2]:<8s}"
                + "".join(f" {dn[j]*1000.0:10.2f}" for j in range(3))
                + "".join(f" {s_comb[j]:10.2f}" for j in range(3))
                + "".join(f" {s_tran[j]:10.2f}" for j in range(3))
                + "\n"
            )

            for j in range(3):
                var_jj = self.sys1_cov[ns1, j, j] + self.sys2_cov[ns2, j, j]
                summ[j] += dn[j] / var_jj
                sumv[j] += dn[j]**2 / var_jj
                sumw[j] += 1.0 / var_jj

        n_fund = len(self.fund_links)
        summ[3] = (summ[0] + summ[1]) / 2.0
        sumv[3] = (sumv[0] + sumv[1]) / 2.0
        sumw[3] = (sumw[0] + sumw[1]) / 2.0

        for j in range(4):
            wmean = summ[j] / sumw[j] * 1000.0 if sumw[j] > 0 else 0.0
            nrms  = np.sqrt(sumv[j] / n_fund) if n_fund > 0 else 0.0
            wrms  = (np.sqrt(1.0 / sumw[j] * n_fund) * nrms * 1000.0
                     if sumw[j] > 0 else 0.0)
            msg = (f"S Component {self._COMP_LABS[j]} # {n_fund:5d}"
                   f" WMean {wmean:6.2f} WRMS {wrms:6.2f}"
                   f" mm/yr, NRMS {nrms:7.3f}")
            print(msg)

    # -----------------------------------------------------------------------
    # Apply transformation to sys1
    # -----------------------------------------------------------------------

    def _update_tran(self) -> None:
        """
        Apply estimated transformation to sys1 velocities and update
        the covariance to include transformation uncertainty.
        Mirrors velrot.f subroutine update_tran.

        Note: the sigma update uses the (intentional) index convention from
        the original Fortran, where cov_neu[0,0] (North transform var) is
        added to sys1_cov[1,1] (East slot) and vice versa.  This matches
        the Fortran behaviour exactly.
        """
        for i in range(len(self.sys1_names)):
            _, loc, R = rotate_geod(self.sys1_coord[i, :, 1], "XYZ", "NEU",
                                    self.sys1_coord[i, :, 0])
            A_neu, A_xyz = get_parts(self.sys1_coord[i, :, 0], R, self.num_parn)

            # Apply transformation in XYZ
            dx = np.zeros(3)
            for k in range(self.num_parn):
                dx += A_xyz[:, k] * self.trans_parm[k]
            self.sys1_coord[i, :, 1] += dx

            # Propagate transformation uncertainty
            cov_neu = var_comp(A_neu, self.norm_eq)  # (3,3) in (m/yr)²

            # Update covariance (Fortran convention: N↔E swap in transform term)
            # dsig(2) = sqrt(sys1_cov(2,2,i) + cov_neu(1,1)) * 1000
            # dsig(1) = sqrt(sys1_cov(1,1,i) + cov_neu(2,2)) * 1000
            # where index 1=N, 2=E in Fortran → index 0=N, 1=E in Python
            dsig_e = np.sqrt(max(0, self.sys1_cov[i, 1, 1] + cov_neu[0, 0])) * 1000.0
            dsig_n = np.sqrt(max(0, self.sys1_cov[i, 0, 0] + cov_neu[1, 1])) * 1000.0
            dsig_u = np.sqrt(max(0, self.sys1_cov[i, 2, 2] + cov_neu[2, 2])) * 1000.0

            rho_ne = ((self.sys1_cov[i, 0, 1] + cov_neu[0, 1]) /
                      (dsig_n * dsig_e / 1.0e6))

            self.sys1_cov[i, 1, 1] = dsig_e**2 / 1.0e6   # E variance
            self.sys1_cov[i, 0, 0] = dsig_n**2 / 1.0e6   # N variance
            self.sys1_cov[i, 2, 2] = dsig_u**2 / 1.0e6   # U variance
            self.sys1_cov[i, 0, 1] = rho_ne * np.sqrt(
                self.sys1_cov[i, 0, 0] * self.sys1_cov[i, 1, 1])
            self.sys1_cov[i, 1, 0] = self.sys1_cov[i, 0, 1]

    # -----------------------------------------------------------------------
    # Output transformed velocity field
    # -----------------------------------------------------------------------

    def _check_cp(self, ns1: int, ns2: int) -> str:
        """
        Return '*' if sys1 site ns1 is within cp_dist of any sys2 site,
        '+' if sys2 site ns2 is within cp_dist of any sys1 site,
        ' ' otherwise.
        Mirrors velrot.f subroutine check_cp.
        """
        if self.cp_dist <= 0.0:
            return " "
        if ns1 >= 0:
            for k in range(len(self.sys2_names)):
                dl = np.linalg.norm(
                    self.sys1_coord[ns1, :, 0] - self.sys2_coord[k, :, 0])
                if dl < self.cp_dist:
                    return "*"
        elif ns2 >= 0:
            for k in range(len(self.sys1_names)):
                dl = np.linalg.norm(
                    self.sys2_coord[ns2, :, 0] - self.sys1_coord[k, :, 0])
                if dl < self.cp_dist:
                    return "+"
        return " "

    def _output_frame(self, fh) -> None:
        """
        Write the transformed sys1 and the sys2-only sites.
        Mirrors velrot.f subroutine output_frame.

        Format 420 (sys1):
            2(1x,f10.5),1x,6(1x,f7.2),1x,f6.3,2x,3(1x,f7.2), 1x,a8,a1
        Format 520 (sys2):
            a,f10.5,1x,F10.5,1x,6(1x,f7.2),1x,f6.3,2x,3(1x,f7.2), 1x,a8,a1
        """
        fh.write(
            "\n* SYSTEM 1 Velocities transformed to SYSTEM 2 \n"
            "*   Long.       Lat.         E & N Rate      E & N Adj."
            "      E & N +-   RHO        H Rate   H adj.    +-  SITE\n"
            "*  (deg)      (deg)           (mm/yr)       (mm/yr)"
            "       (mm/yr)                 (mm/yr)\n"
        )

        for i in range(len(self.sys1_names)):
            dn, loc_coord, _ = rotate_geod(self.sys1_coord[i, :, 1], "XYZ", "NEU",
                                            self.sys1_coord[i, :, 0])
            dn_mm = dn * 1000.0

            lon = np.degrees(loc_coord[1])
            lat = 90.0 - np.degrees(loc_coord[0])

            se  = np.sqrt(max(0, self.sys1_cov[i, 1, 1])) * 1000.0
            sn  = np.sqrt(max(0, self.sys1_cov[i, 0, 0])) * 1000.0
            su  = np.sqrt(max(0, self.sys1_cov[i, 2, 2])) * 1000.0
            rho = (self.sys1_cov[i, 0, 1] /
                   (sn * se / 1.0e6)) if (sn > 0 and se > 0) else 0.0

            sym = self._check_cp(i, -1)

            fh.write(
                f" {lon:10.5f} {lat:10.5f}"
                f" {dn_mm[1]:7.2f} {dn_mm[0]:7.2f}"
                f" {dn_mm[1]:7.2f} {dn_mm[0]:7.2f}"
                f" {se:7.2f} {sn:7.2f}"
                f" {rho:6.3f}  {dn_mm[2]:7.2f} {dn_mm[2]:7.2f}"
                f" {su:7.2f} {self.sys1_names[i]:<8s}{sym}\n"
            )

        # sys2 sites
        fh.write(
            "\n* SYSTEM 2 Velocities except those in SYSTEM 1 \n"
            "*   Long.       Lat.         E & N Rate      E & N Adj."
            "      E & N +-   RHO        H Rate   H adj.    +-  SITE\n"
            "*  (deg)      (deg)           (mm/yr)       (mm/yr)"
            "       (mm/yr)                 (mm/yr)\n"
        )

        sys1_nameset = set(self.sys1_names)

        for i in range(len(self.sys2_names)):
            dn, loc_coord, _ = rotate_geod(self.sys2_coord[i, :, 1], "XYZ", "NEU",
                                            self.sys2_coord[i, :, 0])
            dn_mm = dn * 1000.0

            lon = np.degrees(loc_coord[1])
            lat = 90.0 - np.degrees(loc_coord[0])

            se  = np.sqrt(max(0, self.sys2_cov[i, 1, 1])) * 1000.0
            sn  = np.sqrt(max(0, self.sys2_cov[i, 0, 0])) * 1000.0
            su  = np.sqrt(max(0, self.sys2_cov[i, 2, 2])) * 1000.0
            rho = (self.sys2_cov[i, 0, 1] /
                   (sn * se / 1.0e6)) if (sn > 0 and se > 0) else 0.0

            sym = self._check_cp(-1, i)

            # Leading character: '-' if site is also in sys1, ' ' otherwise
            leader = "-" if self.sys2_names[i] in sys1_nameset else " "

            fh.write(
                f"{leader}{lon:10.5f} {lat:10.5f}"
                f" {dn_mm[1]:7.2f} {dn_mm[0]:7.2f}"
                f" {dn_mm[1]:7.2f} {dn_mm[0]:7.2f}"
                f" {se:7.2f} {sn:7.2f}"
                f" {rho:6.3f}  {dn_mm[2]:7.2f} {dn_mm[2]:7.2f}"
                f" {su:7.2f} {self.sys2_names[i]:<8s}{sym}\n"
            )

    # -----------------------------------------------------------------------
    # Convenience: get results as dicts
    # -----------------------------------------------------------------------

    def get_parameters(self) -> Dict:
        """
        Return transformation parameters in display units.

        Returns
        -------
        dict with keys: Tx, Ty, Tz (mm/yr), Rx, Ry, Rz (mas/yr),
                        Scale (ppb/yr), and sigma_* equivalents.
        """
        labs = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "Scale"]
        return {
            **{labs[i]: self.trans_parm_out[i] for i in range(7)},
            **{f"sigma_{labs[i]}": self.trans_sigma[i] for i in range(7)},
            "rms_mm_yr": self.rms_fit,
            "chi_fit": self.chi_fit,
            "num_data": self.num_data,
            "num_fund": len(self.fund_links),
        }

    def get_sys1_velocities(self) -> List[Dict]:
        """
        Return the transformed System 1 velocities as a list of dicts
        with keys: name, lon, lat, ve, vn, vu, se, sn, su, rho_ne.
        """
        result = []
        for i in range(len(self.sys1_names)):
            dn, loc, _ = rotate_geod(self.sys1_coord[i, :, 1], "XYZ", "NEU",
                                     self.sys1_coord[i, :, 0])
            lon = np.degrees(loc[1])
            lat = 90.0 - np.degrees(loc[0])
            se  = np.sqrt(max(0, self.sys1_cov[i, 1, 1])) * 1000.0
            sn  = np.sqrt(max(0, self.sys1_cov[i, 0, 0])) * 1000.0
            su  = np.sqrt(max(0, self.sys1_cov[i, 2, 2])) * 1000.0
            rho = (self.sys1_cov[i, 0, 1] / (sn * se / 1.0e6)
                   if (sn > 0 and se > 0) else 0.0)
            result.append(dict(name=self.sys1_names[i], lon=lon, lat=lat,
                               ve=dn[1]*1000, vn=dn[0]*1000, vu=dn[2]*1000,
                               se=se, sn=sn, su=su, rho_ne=rho))
        return result


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def main():
    """
    Command-line interface matching velrot runstring:
        python pyvelrot.py sys1 frame1 sys2 frame2 outname out_frame \
                           fund_file height_wght param_opt
    """
    args = sys.argv[1:]
    if len(args) < 5:
        print("Usage: pyvelrot.py sys1 frame1 sys2 frame2 outname "
              "[out_frame] [fund_file] [height_weight] [param_opt]")
        sys.exit(1)

    sys1_file  = args[0]
    sys1_frame = args[1]
    sys2_file  = args[2]
    sys2_frame = args[3]
    out_file   = args[4]
    out_frame  = args[5] if len(args) > 5 else sys1_frame
    fund_file  = args[6] if len(args) > 6 else ""
    h_weight   = float(args[7]) if len(args) > 7 else 1.0
    p_opt      = args[8] if len(args) > 8 else "TR"

    vr = PyVelrot()
    vr.run(sys1_file, sys1_frame, sys2_file, sys2_frame,
           out_file, out_frame, fund_file, h_weight, p_opt)


if __name__ == "__main__":
    main()
