"""Functions used to atomtype a gmso.Topology."""
from operator import attrgetter

import networkx as nx
from foyer.atomtyper import AtomTypingRulesProvider, find_atomtypes
from foyer.smarts import SMARTS
from foyer.topology_graph import TopologyGraph

from gmso import Atom
from gmso.atomtyping.isomorph import (  # check refs
    partition_isomorphic_topology_graphs,
)

__all__ = ["apply"]


def apply(
    top,
    forcefields,
    identify_connected_components=True,
    use_residue_info=False,
    assert_bond_params=True,
    assert_angle_params=True,
    assert_dihedral_params=True,
    assert_improper_params=False,
):
    """Set Topology parameter types from GMSO Forcefields.

    Parameters
    ----------
    top: gmso.core.topology.Topology, required
        The GMSO topology on which to apply forcefields
    forcefields: dict, required
        The keys are labels that match the subtopology label or site residue_name, and the
        values are gmso Forcefield objects that gets applied to the specified molecule
    identify_connected_components: bool, optional, default=True
        A flag to determine whether or not to search the topology for repeated disconnected
        structures, otherwise known as molecules and type each molecule only once.
    use_residue_info: bool, optional, default=False
        A flag to determine whether or not to look at site.residue_name to look parameterize
        each molecule only once. Will only be used if identify_connected_components=False
    assert_bond_params : bool, optional, default=True
        If True, Foyer will exit if parameters are not found for all system
        bonds.
    assert_angle_params : bool, optional, default=True
        If True, Foyer will exit if parameters are not found for all system
        angles.
    assert_dihedral_params : bool, optional, default=True
        If True, Foyer will exit if parameters are not found for all system
        proper dihedrals.
    assert_improper_params : bool, optional, default=False
        If True, Foyer will exit if parameters are not found for all system
        improper dihedrals.
    """
    top.identify_connections()
    top_graph = get_modified_topology_graph(top)
    if identify_connected_components:  # True, (True,False)
        # populate molecule information
        isomorphic_top_graphs = partition_isomorphic_topology_graphs(top_graph)
        apply_atomtypes_from_isomorphic(
            top, isomorphic_top_graphs, forcefields
        )  # what if only one forcefield?, traversing done in this function
    elif not use_residue_info:  # False, False
        # get entire topgraph
        topgraph = get_modified_topology_graph(top)
        # get rules provider from gmso_ff
        # Assumes forcefields is just one forcefield
        ff_of_interest = next(iter(forcefields.values()))
        at_rules_provider = get_atomtyping_rules_provider(gmso_ff=ff_of_interest)
        # create a typemap of that topgraph
        typemap = find_atomtypes(topgraph, at_rules_provider)
        # apply typemap
        traverse_typemap(top, typemap, ff_of_interest)
    elif use_residue_map:  # False, True
        """Not Implemented
        maps = {}
        for key, val in forcefields.items():
            # generate topgraph of items with "key" residue name
            # NOTE CAL: this iter_sites_by method won't work currently, can only iter one at a time.
            res_iter = top.iter_sites_by(["residue_name", "residue_number"], [key, 1]) #Need to note that resnumber index 1
            # create a top/graph from this set of sites
            #subgraph = FUNCTION(res_iter)
            # append that typemap to a dict with residue name as key
            #maps[key] = find_atomtypes(subgraph, val)
        # generate a full typemap by iterating through all sites and copying over map info from `maps`
        #typemap = FUNCTION(top, maps)
        # apply typemap
        traverse_typemap(top, typemap, forcefields)
        """
        raise (
            GMSOError,
            "Using site residue matching to match substructures to a given forcefield is not implemented yet",
        )

    if assert_bond_params:
        apply_connection_types(top, forcefields, "bond")
    if assert_angle_params:
        apply_connection_types(top, forcefields, "angle")
    if assert_dihedral_params:
        apply_connection_types(top, forcefields, "dihedral")
    if assert_improper_params:
        apply_connection_types(top, forcefields, "improper")

    return top


def get_modified_topology_graph(gmso_topology):
    """Return a TopologyGraph with relevant attributes from an GMSO topology.

    Parameters
    ----------
    gmso_topology: gmso.Topology
        The GMSO Topology

    Returns
    -------
    TopologyGraph
        The equivalent TopologyGraph of the openFF Topology `openff_topology`
    """
    top_graph = TopologyGraph()
    for atom in gmso_topology.sites:
        if isinstance(atom, Atom):
            if atom.name.startswith("_"):
                top_graph.add_atom(
                    name=atom.name,
                    index=gmso_topology.get_index(atom),
                    atomic_number=None,
                    element=atom.name,
                    subtopology_label=atom.label,
                )

            else:
                top_graph.add_atom(
                    name=atom.name,
                    index=gmso_topology.get_index(atom),
                    atomic_number=atom.element.atomic_number,
                    element=atom.element.symbol,
                    subtopology_label=atom.label,
                )

    for top_bond in gmso_topology.bonds:
        atoms_indices = [
            gmso_topology.get_index(atom)
            for atom in top_bond.connection_members
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
    parser = SMARTS({})
    for atom_type_name, atom_type in gmso_ff.atom_types.items():
        if atom_type.definition:
            atom_type_defs[atom_type_name] = atom_type.definition
        if atom_type.overrides:
            atom_type_overrides[atom_type_name] = atom_type.overrides

    return AtomTypingRulesProvider(
        atom_type_defs, atom_type_overrides, {}, parser
    )


def apply_atomtypes_from_isomorphic(top, isomorphic_top_graphs, forcefields):
    """Set Topology atomtypes from isomorphic structures matching supplied forcefields.

    Parameters
    ----------
    top: gmso.core.topology.Topology, required
        The GMSO topology on which to apply atomtypes
    isomorphic_top_graphs: dict, required
        The dictionary mapping which breaks up molecules into one graph structure, and a list of
        the corresponding subgraphs in the TopologyGraph that are isomorphic
    forcefields: dict, required
        The keys are labels that match the subtopology labels for each isomorphic graph, and the
        values are gmso Forcefield objects that gets applied to the specified molecule
    """
    for top_graph, identical in isomorphic_top_graphs.items():
        sub_top_label = next(
            iter(
                nx.classes.function.get_node_attributes(
                    top_graph, "atom_data"
                ).values()
            )
        ).subtopology_label
        ff_of_interest = forcefields[sub_top_label]
        at_rules_provider = get_atomtyping_rules_provider(
            gmso_ff=ff_of_interest
        )
        sub_type_map = find_atomtypes(top_graph, at_rules_provider)
        traverse_typemap(top, sub_type_map, ff_of_interest)
        for graph, mapping in identical:
            for node in graph:
                mapped = sub_type_map[mapping[node]]["atom_type"]
                top._sites[node].atom_type = mapped.clone()


def traverse_typemap(top, type_map, forcefield):
    """Take a typemap and applies those atomtypes to the Topology.

    Parameters
    ----------
    top: gmso.core.topology.Topology, required
        The GMSO topology on which to apply atomtypes
    typemap: dict, required
        The dictionary of atom_index with iterated matching to get the atomtype of that site
    forcefield: gmso.forcefield.Forcefield, required
        The GMSO forcefield object which has information associated with each unique atomtype
    """
    for index in type_map:
        site_of_interest = top._sites[index]
        site_of_interest.atom_type = forcefield.atom_types[
            type_map[index]["atomtype"]
        ].clone()
        type_map[index]["atom_type"] = site_of_interest.atom_type


def apply_residue_names(top, more):  # TODO
    """Take a topology and applies residue info to all sites.

    Parameters
    ----------
    top: gmso.core.topology.Topology, required
        The GMSO topology on which to apply atomtypes
    more: TBD, required
        The information necessary to map residue name and number for each site
    """
    return


def apply_connection_types(top, forcefields, parameter):
    """Set Topology connection types from matching atomclass information in forcefield.

    Parameters
    ----------
    top: gmso.core.topology.Topology, required
        The GMSO topology on which to apply paramter_types
    forcefields: dict, required
        The keys are labels that match the subtopology labels for each molecule, and the
        values are gmso Forcefield objects that gets applied to the specified molecule
    parameter: str, required
        The connection type to parameterize. Can any of "bond", "angle", "dihedral", "improper"
    """
    connect_members = {
        "bond": [0, 1],
        "angle": [0, 1, 2],
        "dihedral": [0, 1, 2, 3],
        "improper": [0, 1, 2, 3],
    }  # TODO validate improper order
    for connect in getattr(top, "_" + parameter + "s"):
        sub_top_label = connect.connection_members[0].label_
        ff_of_interest = forcefields[sub_top_label]
        type_getter = attrgetter("atom_type.atomclass")
        member_types = [
            type_getter(getattr(connect, "connection_members")[i])
            for i in connect_members[parameter]
        ]
        ctype = ff_of_interest.get_potential(
            group=parameter + "_type", key=member_types
        ).clone()
        setattr(connect, parameter + "_type", ctype)
