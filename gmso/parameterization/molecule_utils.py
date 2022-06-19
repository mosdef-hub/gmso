"""Utilities for application of a particular forcefield to a molecule."""
from gmso.abc.abstract_site import MoleculeType, ResidueType


def _conn_in_molecule(connection, molecule):
    """Check if all the members in a connection belong to a molecule (namedtuple)."""
    return all(
        site.molecule == molecule for site in connection.connection_members
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


def assert_no_boundary_bonds(top, mol_or_res):
    """Given a molecule tag, assert that no bonds exist between its sites and external sites."""
    namedtuple_types = {MoleculeType: "molecule", ResidueType: "residue"}
    key = namedtuple_types[type(mol_or_res)]
    for bond in top.bonds:
        site_pairs = bond.connection_members
        assertion_msg = (
            "Site {} is in the molecule {}, but its bonded partner {} is not."
        )
        molecule_sites = top.iter_sites(key, mol_or_res)
        if site_pairs[0] in top.iter_sites(key, mol_or_res):
            assert site_pairs[1] in molecule_sites, assertion_msg.format(
                site_pairs[0].name, mol_or_res.name, site_pairs[1].name
            )
        elif site_pairs[1] in top.iter_sites(key, mol_or_res):
            assert site_pairs[0] in molecule_sites, assertion_msg.format(
                site_pairs[1].name, mol_or_res.name, site_pairs[0].name
            )
