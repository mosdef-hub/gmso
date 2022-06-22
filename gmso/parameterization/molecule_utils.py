"""Utilities for application of a particular forcefield to a molecule."""
from gmso.abc.abstract_site import MoleculeType, ResidueType


def _conn_in_molecule(connection, molecule, group=False):
    """Check if all the members in a connection belong to a molecule (namedtuple)."""
    attr = "group" if isinstance(molecule, str) else "molecule"
    return all(
        getattr(site, attr) == molecule
        for site in connection.connection_members
    )


def _molecule_connections(top, molecule, attr):
    """Return all the connections belonging to a molecule."""
    return filter(
        lambda conn: _conn_in_molecule(conn, molecule),
        getattr(top, attr),
    )


def molecule_bonds(top, molecule):
    """Given a molecule (namedtuple), return its bonds."""
    return _molecule_connections(top, molecule, "bonds")


def molecule_angles(top, molecule):
    """Given a molecule (namedtuple), return its angles."""
    return _molecule_connections(top, molecule, "angles")


def molecule_dihedrals(top, molecule):
    """Given a molecule (namedtuple), return its dihedrals."""
    return _molecule_connections(top, molecule, "dihedrals")


def molecule_impropers(top, molecule):
    """Given a molecule (namedtuple), return its impropers."""
    return _molecule_connections(top, molecule, "impropers")


def assert_no_boundary_bonds(top):
    """Assert that all bonds in the topology belongs to only one molecule."""
    assertion_msg = "Site {} is in the molecule {}, but its bonded partner {} is in the molecule {}."
    for bond in top.bonds:
        site1, site2 = bond.connection_members
        assert site1.molecule == site2.molecule, assertion_msg.format(
            site1.name, site1.molecule, site2.name, site2.molecule
        )
