"""Utilities for atomtyping a gmso topology with foyer."""
from collections import namedtuple

from foyer.atomtyper import AtomTypingRulesProvider, find_atomtypes
from foyer.topology_graph import TopologyGraph

from gmso.core.atom import Atom
from gmso.parameterization.subtopology_utils import subtop_bonds
from gmso.utils.decorators import experimental


@experimental
def get_topology_graph(gmso_topology, atomdata_populator=None):
    """Return a TopologyGraph with relevant attributes from an GMSO topology.

    Parameters
    ----------
    gmso_topology: gmso.Topology-like
        The GMSO Topology

    atomdata_populator: callable, default=None
        A function that will be called with the following arguments `gmso_topology` as well as `atom` to pass extra
        arguments to the foyer.topology_graph.AtomData object

    Notes
    -----
    The gmso topology here is a duck type.

    Returns
    -------
    foyer.topology_graph.TopologyGraph
        A light networkx representation of the topology
    """
    top_graph = TopologyGraph()
    for atom in gmso_topology.sites:
        if isinstance(atom, Atom):
            kwargs = (
                atomdata_populator(gmso_topology, atom)
                if atomdata_populator
                else {}
            )
            if atom.name.startswith("_"):
                top_graph.add_atom(
                    name=atom.name,
                    index=gmso_topology.get_index(atom),
                    atomic_number=None,
                    element=atom.name,
                    **kwargs,
                )

            else:
                top_graph.add_atom(
                    name=atom.name,
                    index=gmso_topology.get_index(atom),
                    atomic_number=atom.element.atomic_number,
                    element=atom.element.symbol,
                    **kwargs,
                )

    for top_bond in gmso_topology.bonds:
        atoms_indices = [
            gmso_topology.get_index(atom)
            for atom in top_bond.connection_members
        ]
        top_graph.add_bond(atoms_indices[0], atoms_indices[1])

    return top_graph


@experimental
def get_topology_graph_from_subtop(subtopology):
    """Get an equivalent topology graph for a sub-topology."""
    subtop_named_tuple = namedtuple("subtopology", ("sites", "bonds"))
    return get_topology_graph(
        subtop_named_tuple(subtopology.sites, subtop_bonds(subtopology))
    )


@experimental
def get_atomtyping_rules_provider(gmso_ff):
    """Return a foyer AtomTypingRulesProvider from a GMSO forcefield.

    Parameters
    ----------
    gmso_ff: gmso.core.forcefield.Forcefield
        The GMSO forcefield object to extract the rules from

    Returns
    -------
    AtomTypingRulesProvider
        The foyer.atomtyper.AtomTypingRulesProvider object used to parse atomtype definitions.
        Typically, SMARTS is the ruleset of choice. See https://github.com/mosdef-hub/foyer/issues/63
        for curently supported features in Foyer.
    """
    atom_type_defs = {}
    atom_type_overrides = {}
    for atom_type_name, atom_type in gmso_ff.atom_types.items():
        if atom_type.definition:
            atom_type_defs[atom_type_name] = atom_type.definition
        if atom_type.overrides:
            atom_type_overrides[atom_type_name] = atom_type.overrides

    return AtomTypingRulesProvider(
        atom_type_defs, atom_type_overrides, gmso_ff.non_element_types
    )


@experimental
def typemap_dict(topology_graph, atomtyping_rules_provider, max_iter=10):
    """Return a dictionary of typemap, by finding atomtypes in foyer."""
    return find_atomtypes(topology_graph, atomtyping_rules_provider, max_iter)
