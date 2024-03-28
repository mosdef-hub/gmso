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


def sort_connection_members(connection, sort_by="name", top=None):
    """Sort connection_members of connection.

    Parameters
    ----------
    connection : gmso.Bond, gmso.Angle, gmso.Dihedral, gmso.Improper
        The connection made up of sites to sort by
    sort_by : str, default="name"
        The attribute of the site to sort by. Can take "name", "atom_type",
        and "atom_class" as the sorting attribute.
    top : gmso.Topology
        The topology associated with the connection. This is used in sort_by='index'
        to provide the index of the sites.
    """
    if sort_by == "name":
        sorting_key = lambda site: natural_sort(site.name)
    elif sort_by == "atom_type":
        sorting_key = lambda site: natural_sort(site.atom_type.name)
    elif sort_by == "atomclass":
        sorting_key = lambda site: natural_sort(site.atom_type.atomclass)
    elif sort_by == "index":
        if top is None or connection.connection_members[0] not in top.sites:
            raise ValueError(
                f"Must provide topology associated with {connection}. Provided {top}"
            )
        sorting_key = lambda site: top.get_index(site)
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
    """Get tuple of classes for a topology potential based on member_classes.

    Useful for sorting a long list of potentials by the potential.member_classes.
    Returns the sorted potential based on the ordering of its
    `site.atom_type.atomclass`. For instance, with angles, the central atom is
    alwasys listed second. For dihedrals, the middle two atoms are listed at index
    1 and 2, and the outer two atoms are placed at index 0 (which is bonded to index
    1) and index 3 (which is bonded to atom 4). Finally, impropers are organized
    with the central atom at index0, followed by any combination of the other three
    sites.
    Use `sorted(dihedralTypesList, key=sort_by_classes)` for sorting functionality.

    Parameters
    ----------
    potential : gmso.core.ParametricPotential
        Sort one of the potentials, such as an AtomType, BondType, AngleType,
        DihedralType, or ImproperType

    Returns
    -------
    tuple : sorted potential.member_classes based on ordering of sites in GMSO
        for that particular connection. i.e. impropers specify that the central
        atom of the improper is listed first.
    """
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
            *sorted(potential.member_classes[1:3]),
            potential.member_classes[3],
        )
    return ValueError(
        f"Potential {potential} not one of {potential_attribute_map.values()}"
    )


def sort_by_types(potential):
    """Get tuple of types for a topology potential based on member_types.

    Useful for sorting a long list of potentials by the potential.member_types.
    Returns the sorted potential based on the ordering of its
    `site.atom_type.name`. For instance, with angles, the central atom is
    alwasys listed second. For dihedrals, the middle two atoms are listed at index
    1 and 2, and the outer two atoms are placed at index 0 (which is bonded to index
    1) and index 3 (which is bonded to atom 4). Finally, impropers are organized
    with the central atom at index0, followed by any combination of the other three
    sites.
    Use `sorted(dihedralTypesList, key=sort_by_types)` for sorting functionality.

    Parameters
    ----------
    potential : gmso.core.ParametricPotential
        Sort one of the potentials, such as an AtomType, BondType, AngleType,
        DihedralType, or ImproperType

    Returns
    -------
    tuple : sorted potential.member_types based on ordering of sites in GMSO
        for that particular connection. i.e. impropers specify that the central
        atom of the improper is listed first.
    """
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
            *sorted(potential.member_types[1:3]),
            potential.member_types[3],
        )
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
            [namesList[0], *sorted(namesList[1:])],
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


# reindex molecule
def reindex_molecules(top):
    """Take a Topology and reset the molecule numbers to index from 0."""
    unique_moleculesDict = {}
    for site in top.sites:
        molecule = site.molecule
        if molecule.name in unique_moleculesDict:
            unique_moleculesDict[molecule.name].add(molecule.number)
        else:
            unique_moleculesDict[molecule.name] = {molecule.number}

    offsetDict = {}
    for molecule in unique_moleculesDict:
        min_val = min(unique_moleculesDict[molecule])
        offsetDict[molecule] = min_val

    for site in top.sites:
        mol_num = site.molecule.number
        site.molecule = site.molecule._replace(
            number=mol_num - offsetDict[site.molecule.name]
        )
