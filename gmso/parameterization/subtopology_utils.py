"""Utilities for application of a particular forcefield to a subtopology."""


def _members_in_subtop(connection, subtop):
    """Check if all the members in a connection belong to a subtopology."""
    return all(site in subtop.sites for site in connection.connection_members)


def _subtop_connections(subtop, attr):
    """Return all the connections belonging to a subtopology."""
    return filter(
        lambda conn: _members_in_subtop(conn, subtop),
        getattr(subtop._parent, attr),
    )


def subtop_bonds(subtop):
    """Given a subtopology, return its bonds."""
    return _subtop_connections(subtop, "bonds")


def subtop_angles(subtop):
    """Given a subtopology, return its angles."""
    return _subtop_connections(subtop, "angles")


def subtop_dihedrals(subtop):
    """Given a subtopology, return its dihedrals."""
    return _subtop_connections(subtop, "dihedrals")


def subtop_impropers(subtop):
    """Given a subtopology, return its impropers."""
    return _subtop_connections(subtop, "impropers")


def assert_no_boundary_bonds(subtop):
    """Given a subtopology, assert that no bonds exist between its sites and external sites."""
    for bond in subtop._parent.bonds:
        site_pairs = bond.connection_members
        assertion_msg = "Site {} is in the subtopology {}, but its bonded partner {} is not."

        if site_pairs[0] in subtop.sites:
            assert site_pairs[1] in subtop.sites, assertion_msg.format(
                site_pairs[0].name, subtop.name, site_pairs[1].name
            )
        elif site_pairs[1] in subtop.sites:
            assert site_pairs[0] in subtop.sites, assertion_msg.format(
                site_pairs[1].name, subtop.name, site_pairs[0].name
            )
