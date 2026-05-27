"""
frame_registry.py — GAMIT/GLOBK frame name registry
=====================================================

Python implementation of the GAMIT Fortran subroutine ``frame_to_frame``
(source: ``kf/utils/frame_names/frame_to_fra.f``, last updated TAH 180117).

All rotation-vector components are stored in **deg/Myr** relative to the
NNR-NUVEL-1A reference field.  All values are transcribed directly from the
Fortran ``frame_data`` array; they are not computed or derived.

Usage
-----
    from frame_registry import frame_to_frame, list_frames

    # Rotation vector (rad/yr) that transforms velocities from EURA_I14 to ITRF14
    rot = frame_to_frame("EURA_I14", "ITRF14")

    # Print all registered frames
    list_frames()

API
---
frame_to_frame(sys_frame, out_frame) → np.ndarray[3]
    Return the rotation vector **rot_vec** (rad/yr) such that

        V_out = V_sys + rot_vec × X

    where X is the ECEF position vector (m).

    *sys_frame* and *out_frame* are case-insensitive.
    Passing "NONE" for either argument returns a zero vector.
    Passing an unrecognised name raises ``ValueError`` listing valid names.

list_frames() → None
    Print a table of all registered frames (name, wx, wy, wz in deg/Myr).

FRAME_REGISTRY : dict[str, tuple[float, float, float]]
    Direct access to the registry mapping ``name → (wx, wy, wz)`` deg/Myr.

Notes
-----
*  The Fortran fallback that searches ``./frames.dat`` or
   ``~/gg/tables/frames.dat`` for additional frame definitions is
   **not implemented**.  The 78 built-in frames cover all standard
   GAMIT frame names.
*  The NUVEL-1 scale factor (``:O`` suffix) is not supported.
   All built-in ITRF-referenced frames already use NUVEL-1A scaling.
*  The NNR-NUVEL-1A frame itself is registered as ``NUV-NNR``
   (rotation vector = 0), which is the common reference from which
   all other vectors are measured.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Internal constants
# ---------------------------------------------------------------------------

_PI = 3.1415926535897932


# ---------------------------------------------------------------------------
# Frame registry
#
# All values in deg/Myr, relative to NNR-NUVEL-1A.
# Order and values match the Fortran frame_data / frame_names arrays
# in frame_to_fra.f (GAMIT version ~10.71, TAH 180117 update).
#
# Sources
# -------
# NUVEL-1A plates (PCFC .. SCOTIA):
#   DeMets, C., R.G. Gordon, D.F. Argus and S. Stein (1990),
#   Current plate motions, Geophys. J. Int., 101, 425-478.
# ITRF2000 plates (_I00 suffix):
#   Altamimi, Z., P. Sillard, and C. Boucher (2002),
#   J. Geophys. Res., 107(B10), 2214, doi:10.1029/2001JB000561.
# ITRF2005 plates (_I05 suffix):
#   Altamimi, Z., X. Collilieux, J. Legrand, B. Garayt and C. Boucher (2007),
#   J. Geophys. Res., 112, B09401, doi:10.1029/2007JB004949.
# ITRF2008 plates (_I08 suffix):
#   Altamimi, Z., L. Metivier, and X. Collilieux (2012),
#   J. Geophys. Res., 117, B07402, doi:10.1029/2011JB008930.
# ITRF2014 plates (_I14 suffix):
#   Altamimi, Z., L. Metivier, P. Rebischung, H. Rouby, X. Collilieux (2017),
#   Geophys. J. Int., 209(3), 1906-1912, doi:10.1093/gji/ggx136.
# ARAB_MCC:
#   McClusky, S., R. Reilinger, S. Mahmoud, D. Ben Sari, A. Tealeb (2003),
#   Geophys. J. Int.
# ARAB_M06:
#   Reilinger, R. et al. (2006), J. Geophys. Res., 111(B5), B05411.
# ---------------------------------------------------------------------------

FRAME_REGISTRY = {
    # ---- NUVEL-1A plates (relative to NNR-NUVEL-1A) ----
    "PCFC":     (-0.08652,   0.27731,  -0.57124),
    "COCO":     (-0.59731,  -1.23788,   0.62596),
    "NAZC":     (-0.08778,  -0.49143,   0.55056),
    "CARB":     (-0.01020,  -0.19395,   0.09058),
    "SAFD":     (-0.05947,  -0.08680,  -0.04985),
    "ANTA":     (-0.04704,  -0.09746,   0.21234),
    "INDI":     ( 0.38216,   0.00229,   0.38904),
    "AUST":     ( 0.44914,   0.29358,   0.35993),
    "AFRC":     ( 0.05105,  -0.17756,   0.22471),
    "ARAB":     ( 0.38302,  -0.02985,   0.38732),
    "EURA":     (-0.05621,  -0.13722,   0.18065),
    "NAFD":     ( 0.01478,  -0.20621,  -0.00877),
    "JUAN":     ( 0.29794,   0.49332,  -0.33346),
    "PHIL":     ( 0.57811,  -0.41024,  -0.55405),
    "RIVERA":   (-0.53801,  -1.77388,   0.69041),
    "SCOTIA":   (-0.02349,  -0.15241,  -0.07277),
    # ---- Special / legacy ----
    "NUV-NNR":  ( 0.0,       0.0,       0.0    ),  # NNR-NUVEL-1A (reference)
    "AM-02":    (-0.0178,    0.0128,    0.0049 ),
    "ITRF93":   ( 0.0333,    0.0806,   -0.0056 ),
    "ITRF94":   ( 0.0,       0.0,       0.0    ),
    "GG_PCFC":  (-0.106303,  0.27194,  -0.60755),
    # ---- ITRF2000 plates ----
    "ANTA_I00": (-0.06344,  -0.08870,   0.20364),
    "AUST_I00": ( 0.40071,   0.32958,   0.32834),
    "EURA_I00": (-0.02246,  -0.13607,   0.22041),
    "NOAM_I00": ( 0.02307,  -0.19187,  -0.01703),
    "PCFC_I00": (-0.10015,   0.27228,  -0.59949),
    "SOAM_I00": (-0.07388,  -0.07484,  -0.04134),
    "ITRF00":   ( 0.0,       0.0,       0.0    ),
    # ---- ARAB special poles ----
    "ARAB_MCC": ( 0.31451,  -0.02398,   0.40449),  # McClusky et al. 2003
    "ARAB_M06": ( 0.33478,  -0.01723,   0.42398),  # Reilinger et al. 2006
    # ---- ITRF2005 plates ----
    "AMUR_I05": (-0.0331,   -0.1457,    0.2237 ),
    "ANTA_I05": (-0.0648,   -0.0915,    0.1928 ),
    "ARAB_I05": ( 0.3735,    0.0331,    0.4412 ),
    "AUST_I05": ( 0.4214,    0.3218,    0.3366 ),
    "CARB_I05": (-0.0460,   -0.1807,    0.1527 ),
    "EURA_I05": (-0.0151,   -0.1439,    0.2172 ),
    "INDI_I05": ( 0.3677,    0.1474,    0.4691 ),
    "NAZC_I05": (-0.0899,   -0.4442,    0.4548 ),
    "NOAM_I05": ( 0.0087,   -0.1913,   -0.0144 ),
    "NUBI_I05": ( 0.0226,   -0.1716,    0.2059 ),
    "OKHT_I05": (-0.0479,   -0.0515,   -0.0440 ),
    "PCFC_I05": (-0.1221,    0.2895,   -0.6053 ),
    "SOAM_I05": (-0.0739,   -0.0892,   -0.0350 ),
    "SOMA_I05": ( 0.0015,   -0.1831,    0.2489 ),
    "YANG_I05": (-0.0533,   -0.1484,    0.2669 ),
    "ITRF05":   ( 0.0,       0.0,       0.0    ),
    # ---- ITRF2008 (Altamimi et al. 2012) ----
    "ITRF08":   ( 0.0,       0.0,       0.0    ),
    "NUBI_I08": ( 0.02639,  -0.16610,   0.20080),
    "ANTA_I08": (-0.07000,  -0.08389,   0.17860),
    "ARAB_I08": ( 0.33390,  -0.01500,   0.41250),
    "AUST_I08": ( 0.41780,   0.32560,   0.34110),
    "CARB_I08": ( 0.01361,  -0.30220,   0.18440),
    "EURA_I08": (-0.02306,  -0.14830,   0.20830),
    "INDI_I08": ( 0.34220,   0.08417,   0.42780),
    "NAZC_I08": (-0.09167,  -0.43080,   0.45140),
    "NOAM_I08": ( 0.009722, -0.18390,  -0.027780),
    "PCFC_I08": (-0.11420,   0.28780,  -0.60170),
    "SOAM_I08": (-0.06750,  -0.08639,  -0.042780),
    "SUND_I08": ( 0.01306,  -0.27780,   0.27080),
    "SOMA_I08": (-0.022220, -0.20690,   0.24920),
    "AMUR_I08": (-0.052780, -0.12280,   0.25420),
    # ---- ITRF2014 (Altamimi et al. 2017) ----
    "ANTA_I14": (-0.0688,   -0.0900,    0.1874 ),
    "ARAB_I14": ( 0.3205,   -0.0378,    0.4011 ),
    "AUST_I14": ( 0.4194,    0.3284,    0.3375 ),
    "EURA_I14": (-0.0235,   -0.1476,    0.2140 ),
    "INDI_I14": ( 0.3205,   -0.0014,    0.4038 ),
    "NAZC_I14": (-0.0925,   -0.4290,    0.4508 ),
    "NOAM_I14": ( 0.0066,   -0.1928,   -0.0176 ),
    "NUBI_I14": ( 0.0274,   -0.1704,    0.2037 ),
    "PCFC_I14": (-0.1135,    0.2907,   -0.6025 ),
    "SOAM_I14": (-0.0751,   -0.0835,   -0.0389 ),
    "SOMA_I14": (-0.0336,   -0.2206,    0.2454 ),
    # ---- Aliases ----
    "ITRF14":   ( 0.0,       0.0,       0.0    ),  # reference for ITRF2014-PMM
    "IGS14":    ( 0.0,       0.0,       0.0    ),  # IGS14 ≡ ITRF14
    "NAM14":    ( 0.0066,   -0.1928,   -0.0176 ),  # NAM14 = NOAM_I14
    "ANT14":    (-0.0688,   -0.0900,    0.1874 ),  # ANT14 = ANTA_I14
    "NAM08":    ( 0.009722, -0.18390,  -0.027780),  # NAM08 = NOAM_I08
    "IGS08":    ( 0.0,       0.0,       0.0    ),  # IGS08 ≡ ITRF08
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def frame_to_frame(sys_frame: str, out_frame: str) -> np.ndarray:
    """
    Return the rotation vector (rad/yr) that maps velocities from *sys_frame*
    to *out_frame*.

    The result satisfies::

        V_out = V_sys + rot_vec × X

    where X is the ECEF site position (m).

    Parameters
    ----------
    sys_frame : str
        Name of the source frame (e.g. ``"EURA_I14"``).
        Case-insensitive.  ``"NONE"`` → zero vector (no rotation).
    out_frame : str
        Name of the target frame.  Same rules as *sys_frame*.

    Returns
    -------
    rot_vec : np.ndarray, shape (3,)
        Rotation vector in rad/yr (X, Y, Z components).

    Raises
    ------
    ValueError
        If *sys_frame* or *out_frame* is not ``"NONE"`` and not found in
        :data:`FRAME_REGISTRY`.  The error message lists all valid names.
    """
    sf = sys_frame.strip().upper()
    of = out_frame.strip().upper()

    # NONE → no-op (matches Fortran behaviour)
    if sf == "NONE" or of == "NONE":
        return np.zeros(3)

    # Same frame → no-op
    if sf == of:
        return np.zeros(3)

    # LIST → print table and return zeros (mirrors Fortran LIST option)
    if sf == "LIST":
        list_frames()
        return np.zeros(3)

    # Look up both frames
    missing = [name for name in (sf, of) if name not in FRAME_REGISTRY]
    if missing:
        valid = ", ".join(sorted(FRAME_REGISTRY))
        raise ValueError(
            f"Unknown frame name(s): {missing}.\n"
            f"Valid names: {valid}\n"
            f"Use list_frames() to print all frames with their rotation vectors."
        )

    ws = np.array(FRAME_REGISTRY[sf], dtype=float)
    wo = np.array(FRAME_REGISTRY[of], dtype=float)

    # Convert deg/Myr → rad/yr  (matches Fortran: * pi / 180.d6)
    return (ws - wo) * (_PI / 180.0 / 1.0e6)


def list_frames() -> None:
    """Print a table of all registered frames (name, wx, wy, wz in deg/Myr)."""
    print(f"* Frame registry: {len(FRAME_REGISTRY)} frames available")
    print(f"* {'#':>3}  {'Name':<10}  {'wx (deg/Myr)':>14}  "
          f"{'wy (deg/Myr)':>14}  {'wz (deg/Myr)':>14}")
    print("* " + "-" * 62)
    for i, (name, (wx, wy, wz)) in enumerate(FRAME_REGISTRY.items(), start=1):
        print(f"* {i:3d}  {name:<10}  {wx:14.6f}  {wy:14.6f}  {wz:14.6f}")
    print("*")
    print("* All values in deg/Myr relative to NNR-NUVEL-1A.")
    print("* Use frame_to_frame(sys, out) to compute rot_vec in rad/yr.")


# ---------------------------------------------------------------------------
# Convenience: available frame names as a frozenset
# ---------------------------------------------------------------------------

FRAME_NAMES = frozenset(FRAME_REGISTRY)
