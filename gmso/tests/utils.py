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


def match_connection_parameters(
    from_gmso_top, from_parmed_top, connection_type="bonds"
):
    """Match connection parameters between two gmso topologies."""
    connection_types_original = {}
    connection_types_mirror = {}
    for connection in getattr(from_parmed_top, connection_type):
        connection_types_mirror[
            tuple(
                from_parmed_top.get_index(member)
                for member in connection.connection_members
            )
        ] = connection

    for connection in getattr(from_gmso_top, connection_type):
        connection_types_original[
            tuple(
                from_gmso_top.get_index(member)
                for member in connection.connection_members
            )
        ] = connection

    for key in connection_types_original:
        conn = connection_types_original[key]
        conn_mirror = connection_types_mirror[key]
        conn_type_attr = connection_type[:-1] + "_type"
        conn_type_mirror = getattr(conn_mirror, conn_type_attr)
        conn_type = getattr(conn, conn_type_attr)
        for param in conn_type.parameters:
            assert u.allclose_units(
                conn_type_mirror.parameters[param],
                conn_type.parameters[param],
            )
