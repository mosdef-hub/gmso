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


def moit(points, masses, center=np.zeros(3)):
    """Calculate moment of inertia tensor (moit) for rigid bodies.

    Assumes rigid body center is at origin unless center is provided.
    Only calculates diagonal elements.

    Parameters
    ----------
    points : numpy.ndarray (N,3)
        x, y, and z coordinates of the rigid body constituent particles
    masses : numpy.ndarray (N,)
        Masses of the constituent particles
    center : numpy.ndarray (3,), default np.array([0,0,0])
        x, y, and z coordinates of the rigid body center

    Returns
    -------
    numpy.ndarray (3,)
        moment of inertia tensor for the rigid body center
    """
    points -= center
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    I_xx = np.sum((y**2 + z**2) * masses)
    I_yy = np.sum((x**2 + z**2) * masses)
    I_zz = np.sum((x**2 + y**2) * masses)
    return np.array((I_xx, I_yy, I_zz))
