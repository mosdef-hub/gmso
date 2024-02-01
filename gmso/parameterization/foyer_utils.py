"""Utilities for atomtyping a gmso topology with foyer."""

from collections import namedtuple

from foyer.atomtyper import AtomTypingRulesProvider, find_atomtypes
from foyer.exceptions import FoyerError
from foyer.topology_graph import TopologyGraph

from gmso.core.atom import Atom
from gmso.parameterization.molecule_utils import molecule_bonds


def get_topology_graph(
    gmso_topology, label_type=None, label=None, atomdata_populator=None
):
    """Return a TopologyGraph with relevant attributes from an GMSO topology.

    Parameters
    ----------
    gmso_topology: gmso.Topology-like
        The GMSO Topology
    label_type: str, optional, default=None
        The type of label used to query the sites-group of interest. Accepted options include "group" and "molecule"
    label: str, optional, default=None
        The label used to query the sites-group of interest
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
    atom_index_map = {}

    if label_type:
        assert label_type in ("group", "molecule"), label_type
        is_group = True if label_type == "group" else False
        pseudo_top = namedtuple("PseudoTop", ("sites", "bonds"))
        gmso_topology = pseudo_top(
            tuple(gmso_topology.iter_sites(label_type, label)),
            tuple(molecule_bonds(gmso_topology, label, is_group)),
        )

    if len(gmso_topology.sites) == 0:
        raise FoyerError(
            "Cannot create a topology graph from a topology with no sites."
        )

    for j, atom in enumerate(gmso_topology.sites):
        atom_index_map[id(atom)] = j
        if isinstance(atom, Atom):
            kwargs = (
                atomdata_populator(gmso_topology, atom)
                if atomdata_populator
                else {}
            )
            if atom.name.startswith("_") or not atom.element:
                top_graph.add_atom(
                    name=atom.name,
                    index=j,  # Assumes order is preserved
                    atomic_number=None,
                    symbol=atom.name,
                    group=atom.group,
                    molecule=atom.molecule.name if atom.molecule else None,
                    **kwargs,
                )
            else:
                top_graph.add_atom(
                    name=atom.name,
                    index=j,  # Assumes order is preserved
                    atomic_number=atom.element.atomic_number,
                    symbol=atom.element.symbol,
                    group=atom.group,
                    molecule=atom.molecule.name if atom.molecule else None,
                    **kwargs,
                )

    for top_bond in gmso_topology.bonds:
        atoms_indices = [
            atom_index_map[id(atom)] for atom in top_bond.connection_members
        ]
        top_graph.add_bond(atoms_indices[0], atoms_indices[1])

    return top_graph


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


def typemap_dict(topology_graph, atomtyping_rules_provider, max_iter=10):
    """Return a dictionary of typemap, by finding atomtypes in foyer."""
    return find_atomtypes(topology_graph, atomtyping_rules_provider, max_iter)
