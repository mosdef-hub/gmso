"""Sorting utilities."""
import re

import gmso
from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType
from gmso.core.dihedral import Dihedral
from gmso.core.dihedral_type import DihedralType
from gmso.core.improper import Improper
from gmso.core.improper_type import ImproperType

potential_attribute_map = {
    Atom: "atom_type",
    Bond: "bond_type",
    Angle: "angle_type",
    Dihedral: "dihedral_type",
    Improper: "improper_type",
}


def _atoi(text):
    """Convert a string to an int."""
    return int(text) if text.isdigit() else text


def natural_sort(text):
    """Given an alphanumeric string, sort using the natural sort algorithm."""
    return [_atoi(a) for a in re.split(r"(\d+)", text)]


def sort_connection_members(connection, sort_by="name"):
    """Sort connection_members of connection."""
    if sort_by == "name":
        sorting_key = lambda site: natural_sort(site.name)
    elif sort_by == "atom_type":
        sorting_key = lambda site: natural_sort(site.atom_type.name)
    elif sort_by == "atomclass":
        sorting_key = lambda site: natural_sort(site.atom_type.atomclass)
    else:
        raise ValueError("Unsupported sort_by value provided.")

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
        if sorting_key(site2) > sorting_key(site3) or (
            sorting_key(site2) == sorting_key(site3)
            and sorting_key(site1) > sorting_key(site4)
        ):
            return [site4, site3, site2, site1]
        else:
            return [site1, site2, site3, site4]
    elif isinstance(connection, gmso.Improper):
        site1, site2, site3, site4 = connection.connection_members
        site2, site3, site4 = sorted([site2, site3, site4], key=sorting_key)
        return [site1, site2, site3, site4]
    else:
        raise TypeError("Provided connection not supported.")


def sort_by_classes(potential):
    """Get list of classes for a topology potential based on memberclass."""
    if isinstance(potential, AtomType):
        return potential.atom_type.atomclass
    elif isinstance(potential, BondType):
        return tuple(sorted(potential.member_classes))
    elif isinstance(potential, AngleType):
        if potential.member_classes[0] > potential.member_classes[2]:
            return tuple(reversed(potential.member_classes))
        else:
            return potential.member_classes
    elif isinstance(potential, DihedralType):
        if potential.member_classes[1] > potential.member_classes[2] or (
            potential.member_classes[1] == potential.member_classes[2]
            and potential.member_classes[0] > potential.member_classes[3]
        ):
            return tuple(reversed(potential.member_classes))
        else:
            return potential.member_classes
    elif isinstance(potential, ImproperType):
        return (
            potential.member_classes[0],
            *potential.member_classes[1:],
        )  # could sort using `sorted`
    return ValueError(
        f"Potential {potential} not one of {potential_attribute_map.values()}"
    )


def sort_by_types(potential):
    """Get list of types for a topology potential based on membertype."""
    if isinstance(potential, AtomType):
        return potential.name
    elif isinstance(potential, BondType):
        return tuple(sorted(potential.member_types))
    elif isinstance(potential, AngleType):
        if potential.member_types[0] > potential.member_types[2]:
            return tuple(reversed(potential.member_types))
        else:
            return potential.member_types
    elif isinstance(potential, DihedralType):
        if potential.member_types[1] > potential.member_types[2] or (
            potential.member_types[1] == potential.member_types[2]
            and potential.member_types[0] > potential.member_types[3]
        ):
            return tuple(reversed(potential.member_types))
        else:
            return potential.member_types
    elif isinstance(potential, ImproperType):
        return (
            potential.member_types[0],
            *potential.member_types[1:],
        )  # could sort using `sorted`
    return ValueError(
        f"Potential {potential} not one of {potential_attribute_map.values()}"
    )


def sort_connection_strings(namesList, improperBool=False):
    """Sort list of strings for a connection to get proper ordering of the connection.

    Parameters
    ----------
    namesList : list
        List of strings connected to a compound to sort.
    improperBool : bool, option,  default=False
        whether or not a four member list refers to an improper
    """
    if len(namesList) == 2:  # assume bonds
        return tuple(sorted(namesList))
    elif len(namesList) == 3:
        if namesList[0] > namesList[2]:
            return tuple(reversed(namesList))
        else:
            return tuple(namesList)
    elif len(namesList) == 4 and improperBool:
        return tuple(
            namesList[0],
            sorted(*namesList[1:]),
        )
    elif len(namesList) == 4 and not improperBool:
        if namesList[1] > namesList[2] or (
            namesList[1] == namesList[2] and namesList[0] > namesList[3]
        ):
            return tuple(reversed(namesList))
        else:
            return tuple(namesList)
    else:
        return ValueError(
            f"Cannot sort {namesList}. It is not a length of 2,3, or 4 members."
        )
