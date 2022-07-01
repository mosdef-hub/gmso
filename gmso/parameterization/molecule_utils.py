"""Utilities for application of a particular forcefield to a molecule."""


def _conn_in_molecule(connection, label, is_group=False):
    """Check if all the members in a connection belong to a molecule (namedtuple)."""
    if is_group:
        return all(
            getattr(site, "group") == label
            for site in connection.connection_members
        )
    else:
        if isinstance(label, str):
            return all(
                getattr(site, "molecule").name == label
                for site in connection.connection_members
            )
        else:
            return all(
                getattr(site, "molecule") == label
                for site in connection.connection_members
            )


def _molecule_connections(top, molecule, attr, is_group=False):
    """Return all the connections belonging to a molecule."""
    return filter(
        lambda conn: _conn_in_molecule(conn, molecule, is_group),
        getattr(top, attr),
    )


def molecule_bonds(top, molecule, is_group=False):
    """Given a molecule (namedtuple), return its bonds."""
    return _molecule_connections(top, molecule, "bonds", is_group)


def molecule_angles(top, molecule, is_group=False):
    """Given a molecule (namedtuple), return its angles."""
    return _molecule_connections(top, molecule, "angles", is_group)


def molecule_dihedrals(top, molecule, is_group=False):
    """Given a molecule (namedtuple), return its dihedrals."""
    return _molecule_connections(top, molecule, "dihedrals", is_group)


def molecule_impropers(top, molecule, is_group=False):
    """Given a molecule (namedtuple), return its impropers."""
    return _molecule_connections(top, molecule, "impropers", is_group)


def assert_no_boundary_bonds(top):
    """Assert that all bonds in the topology belongs to only one molecule."""
    assertion_msg = "Site {} is in the molecule {}, but its bonded partner {} is in the molecule {}."
    for bond in top.bonds:
        site1, site2 = bond.connection_members
        assert site1.molecule == site2.molecule, assertion_msg.format(
            site1.name, site1.molecule, site2.name, site2.molecule
        )
