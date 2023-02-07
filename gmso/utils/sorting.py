"""Sorting utilities."""
import re

import gmso


def _atoi(text):
    """Convert a string to an int."""
    return int(text) if text.isdigit() else text


def natural_sort(text):
    """Given an alphanumeric string, sort using the natural sort algorithm."""
    return [_atoi(a) for a in re.split(r"(\d+)", text)]


def sort_member_types(connection_type):
    """Sort connection_members of connection_type."""
    if isinstance(connection_type, gmso.BondType):
        type1, type2 = connection_type.member_types
        type1, type2 = sorted([type1, type2], key=natural_sort)
        return [type1, type2]
    elif isinstance(connection_type, gmso.AngleType):
        type1, type2, type3 = connection_type.member_types
        type1, type3 = sorted([type1, type3], key=natural_sort)
        return [type1, type2, type3]
    elif isinstance(connection_type, gmso.DihedralType):
        type1, type2, type3, type4 = connection_type.member_types
        if [type2, type3] == sorted([type2, type3], key=natural_sort):
            return [type1, type2, type3, type4]
        else:
            return [type4, type3, type2, type1]
    elif isinstance(connection_type, gmso.ImproperType):
        type1, type2, type3, type4 = connection_type.member_types
        type2, type3, type4 = sorted([type2, type3, type4], key=natural_sort)
        return [type1, type2, type3, type4]
    else:
        raise TypeError("Provided connection_type not supported.")


def sort_connection_members(connection, sort_by="name"):
    """Sort connection_members of connection."""
    if sort_by == "name":

        def sorting_key(site):
            return site.name

    elif sort_by == "atom_type":

        def sorting_key(site):
            return site.atom_type.name

    if isinstance(connection, gmso.Bond):
        site1, site2 = connection.connection_members
        site1, site2 = sorted([site1, site2], key=sorting_key)
        return [site1, site2]
    elif isinstance(connection, gmso.Angle):
        site1, site2, site3 = connection.connection_members
        site1, site3 = sorted([site1, site3], key=sorting_key)
        return [site1, site2, site3]
    elif isinstance(connection, gmso.Dihedral):
        site1, site2, site3, site4 = connection.connection_members
        if [site2, site3] == sorted([site2, site3], key=sorting_key):
            return [site1, site2, site3, site4]
        else:
            return [site4, site3, site2, site1]
    elif isinstance(connection, gmso.Improper):
        site1, site2, site3, site4 = connection.connection_members
        site2, site3, site4 = sorted([site2, site3, site4], key=sorting_key)
        return [site1, site2, site3, site4]
    else:
        raise TypeError("Provided connection not supported.")
