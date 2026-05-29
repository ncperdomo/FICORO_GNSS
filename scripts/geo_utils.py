"""Vectorised geodetic utility functions used across the FICORO_GNSS pipeline."""

import numpy as np


def haversine(lon1, lat1, lon2, lat2):
    """Compute the Haversine great-circle distance (km) between two points.

    All inputs may be scalars or NumPy arrays of the same shape.

    Parameters
    ----------
    lon1, lat1 : float or ndarray  Longitude / latitude of point 1 (degrees).
    lon2, lat2 : float or ndarray  Longitude / latitude of point 2 (degrees).

    Returns
    -------
    float or ndarray
        Distance in kilometres.
    """
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))
