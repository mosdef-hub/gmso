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


def moit(xyz, masses, center=np.zeros(3)):
    """Find the moment of intertia tensor given a set of
    particle coordinates and their corresponding masses.

    This method is used in setting rigid body moments
    of intertia.

    Parameters
    ----------
    xyz : numpy.ndarray (N,3)
        Coordinates of the particles
    masses : numpy.ndarray (N,)
        Masses of the particles
    center : numpy.ndarray (3,), default (0,0,0)
        Coordinates of the particle's center

    Returns
    -------
    numpy.ndarray (3,)
        Moment of inertia tensor for the set of particles.
    """
    xyz -= center
    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]
    Ixx = np.sum((y**2 + z**2) * masses)
    Iyy = np.sum((x**2 + z**2) * masses)
    Izz = np.sum((x**2 + y**2) * masses)
    return np.array((Ixx, Iyy, Izz))
