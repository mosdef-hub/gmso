import os

import unyt as u


def get_path(filename):
    """Given a test filename return its path"""
    _path = os.path.join(os.path.split(__file__)[0], "files", filename)
    return _path


def allclose_units_mixed(u_iter1, u_iter2):
    """Check if array of quantities with mixed dimensions are equivalent.

    Notes
    -----
    The two iterables provided must contain same number of quantities and
    should be able to passed to Python zip function.

    Parameters
    ----------
    u_iter1: list or iterable of u.unit_quantity
        The first iterable/list of unit quantities
    u_iter2: list or iterable of u.unit_quantity
        The second iterable/list of unit quantities

    Returns
    -------
    bool
        True if iter1 is equivalent to iter2
    """
    for q1, q2 in zip(u_iter1, u_iter2):
        if not u.allclose_units(q1, q2):
            return False
    return True
