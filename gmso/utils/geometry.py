"""General module for various geometrical operations."""

import numpy as np


def coord_shift(xyz, box_lengths):
    """Ensure that coordinates are -L/2, L/2.

    Checks if coordinates are -L/2, L/2 and then shifts coordinates
    if necessary. For example, if coordinates are 0, L, then a shift
    is applied to move coordinates to -L/2, L/2. If a shift is not
    necessary, the points are returned unmodified.

    Parameters
    ----------
    xyz : unyt_array of points with shape N x 3
    box : gmso.Box

    Returns
    -------
    xyz : unyt_array of points with shape N x 3
    """
    box_max = box_lengths / 2.0
    box_min = -box_max
    # Shift all atoms
    if np.greater(xyz, box_max).any():
        xyz -= box_max
    elif np.less(xyz, box_min).any():
        xyz += box_max

    return xyz
