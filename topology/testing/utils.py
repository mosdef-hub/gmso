import warnings

import numpy as np


def allclose(a, b, rtol=None, atol=None):
    """Compare two unyt arrays."""
    if a.units != b.units:
        common_unit = _infer_common_unit(a, b)
    else:
        common_unit = a.units

    a = a.in_units(common_unit)
    b = b.in_units(common_unit)

    if atol is None:
        atol = 1e-8 * common_unit

    if rtol is None:
        rtol = 1e-5 * abs(b)

    return (abs(a - b) <= (atol + rtol)).all()


def _infer_common_unit(a, b):
    """Attempt to pick unit best equipped to represent values on two arrays.

    Parameters
    ----------
    a : unyt_array
    b : unyt_array

    Returns
    -------
    common_unit : u.Unit
    """

    if a.units == b.units:
        return a.units

    # Estimate a representative order of magnitude of the values inputs
    order_of_a = np.log10(np.mean(np.abs(a.in_base())))
    order_of_b = np.log10(np.mean(np.abs(b.in_base())))

    if abs(order_of_a - order_of_b) > 6:
        warnings.warn('Values appear to be many orders of magnitude different')

    # Find which array has a typical value nearest 10 ** 0

    if abs(order_of_a) < abs(order_of_b):
        common_unit = a.units
    else:
        common_unit = b.units

    return common_unit
