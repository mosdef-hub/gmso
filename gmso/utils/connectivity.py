"""Module supporting various connectivity methods and operations."""

import networkx as nx
import numpy as np
from boltons.setutils import IndexedSet
from networkx.algorithms import shortest_path_length

from gmso.core.angle import Angle
from gmso.core.dihedral import Dihedral
from gmso.core.improper import Improper

CONNS = {"angle": Angle, "dihedral": Dihedral, "improper": Improper}

EDGES = {
    "angle": ((0, 1),),
    "dihedral": ((0, 1), (1, 2)),
    "improper": ((0, 1), (0, 2), (1, 2)),
}


def identify_connections(top, index_only=False):
    """Identify all possible connections within a topology.

    Parameters
    ----------
    top: gmso.Topology
        The gmso topology for which to identify connections for
    index_only: bool, default=False
        If True, return atom indices that would form the actual connections
        rather than adding the connections to the topology

    Notes: We are using networkx graph matching to match
    the topology's bonding graph to smaller sub-graphs that
    correspond to an angle, dihedral, improper etc.
    The matching is actually done via the line-graph of the
    topology bonding graph.
    [ahy]: IIRC we chose to use line-graph as opposed the actual graph
    because the graph-matching (on the actual graph) would miss certain
    angles/dihedrals/impropers if there were cycles or bridge bonds
    that would effectively hide the angle/dihedral/dihedral
    [ahy]: In the event of virtual sites/drude particles, the matching
    process may have to change in the _detect, _format, or _add methods.
    Personally, I think modifying the _add methods to exclude angles/dihedrals
    with virtual sites would be be the best approach. I
    don't think we would want to change how we construct any of the
    NetworkX graphs.
    """
    compound = nx.Graph()

    for b in top.bonds:
        compound.add_edge(b.connection_members[0], b.connection_members[1])

    compound_line_graph = nx.line_graph(compound)

    angle_matches = _detect_connections(compound_line_graph, top, type_="angle")
    dihedral_matches = _detect_connections(
        compound_line_graph, top, type_="dihedral"
    )
    improper_matches = _detect_connections(
        compound_line_graph, top, type_="improper"
    )

    if not index_only:
        for conn_matches, conn_type in zip(
            (angle_matches, dihedral_matches, improper_matches),
            ("angle", "dihedral", "improper"),
        ):
            if conn_matches:
                _add_connections(top, conn_matches, conn_type=conn_type)
    else:
        return {
            "angles": angle_matches,
            "dihedrals": dihedral_matches,
            "impropers": improper_matches,
        }

    return top


def _add_connections(top, matches, conn_type):
    """Add connections to the topology."""
    for sorted_conn in matches:
        to_add_conn = CONNS[conn_type](
            connection_members=[top.sites[idx] for idx in sorted_conn]
        )
        top.add_connection(to_add_conn, update_types=False)


def _detect_connections(compound_line_graph, top, type_="angle"):
    """Detect available connections in the topology based on bonds."""
    connection = nx.Graph()
    for edge in EDGES[type_]:
        assert len(edge) == 2, "Edges should be of length 2"
        connection.add_edge(edge[0], edge[1])

    matcher = nx.algorithms.isomorphism.GraphMatcher(
        compound_line_graph, connection
    )

    formatter_fns = {
        "angle": _format_subgraph_angle,
        "dihedral": _format_subgraph_dihedral,
        "improper": _format_subgraph_improper,
    }

    conn_matches = IndexedSet()
    for m in matcher.subgraph_isomorphisms_iter():
        new_connection = formatter_fns[type_](m, top)
        conn_matches.add(new_connection)

    if conn_matches:
        conn_matches = _trim_duplicates(conn_matches)

    # Do more sorting of individual connection
    sorted_conn_matches = list()
    for match in conn_matches:
        if type_ in ("angle", "dihedral"):
            if match[0] < match[-1]:
                sorted_conn = match
            else:
                sorted_conn = match[::-1]
        elif type_ == "improper":
            sorted_conn = [match[0]] + sorted(match[1:])
        sorted_conn_matches.append(sorted_conn)

    # Final sorting the whole list
    if type_ == "angle":
        return sorted(
            sorted_conn_matches,
            key=lambda angle: (
                angle[1],
                angle[0],
                angle[2],
            ),
        )
    elif type_ == "dihedral":
        return sorted(
            sorted_conn_matches,
            key=lambda dihedral: (
                dihedral[1],
                dihedral[2],
                dihedral[0],
                dihedral[3],
            ),
        )
    elif type_ == "improper":
        return sorted(
            sorted_conn_matches,
            key=lambda improper: (
                improper[0],
                improper[1],
                improper[2],
                improper[3],
            ),
        )


def _get_sorted_by_n_connections(m):
    """Return sorted by n connections for the matching graph."""
    small = nx.Graph()
    for k, v in m.items():
        small.add_edge(k[0], k[1])
    return sorted(small.adj, key=lambda x: len(small[x])), small


def _format_subgraph_angle(m, top):
    """Format the angle subgraph.

    Since we are matching compound line graphs,
    back out the actual nodes, not just the edges

    Parameters
    ----------
    m : dict
        keys are the compound line graph nodes
        Values are the sub-graph matches (to the angle, dihedral, or improper)
    top : gmso.Topology
        The original Topology

    Returns
    -------
    connection : list of nodes, in order of bonding
        (start, middle, end)
    """
    (sort_by_n_connections, _) = _get_sorted_by_n_connections(m)
    ends = sorted(
        [sort_by_n_connections[0], sort_by_n_connections[1]],
        key=lambda x: top.get_index(x),
    )
    middle = sort_by_n_connections[2]
    return (
        top.get_index(ends[0]),
        top.get_index(middle),
        top.get_index(ends[1]),
    )


def _format_subgraph_dihedral(m, top):
    """Format the dihedral subgraph.

    Since we are matching compound line graphs,
    back out the actual nodes, not just the edges

    Parameters
    ----------
    m : dict
        keys are the compound line graph nodes
        Values are the sub-graph matches (to the angle, dihedral, or improper)
    top : gmso.Topology
        The original Topology

    Returns
    -------
    connection : list of nodes, in order of bonding
        (start, mid1, mid2, end)
    """
    (sort_by_n_connections, small) = _get_sorted_by_n_connections(m)
    start = sort_by_n_connections[0]
    if sort_by_n_connections[2] in small.neighbors(start):
        mid1 = sort_by_n_connections[2]
        mid2 = sort_by_n_connections[3]
    else:
        mid1 = sort_by_n_connections[3]
        mid2 = sort_by_n_connections[2]

    end = sort_by_n_connections[1]
    return (
        top.get_index(start),
        top.get_index(mid1),
        top.get_index(mid2),
        top.get_index(end),
    )


def _format_subgraph_improper(m, top):
    """Format the improper dihedral subgraph.

    Since we are matching compound line graphs,
    back out the actual nodes, not just the edges

    Parameters
    ----------
    m : dict
        keys are the compound line graph nodes
        Values are the sub-graph matches (to the angle, dihedral, or improper)
    top : gmso.Topology
        The original Topology

    Returns
    -------
    connection : list of nodes, in order of bonding
        (central, branch1, branch2, branch3)

    Notes
    -----
    Given the way impropers are matched, sometimes a cyclic 3-ring system gets returned
    """
    (sort_by_n_connections, _) = _get_sorted_by_n_connections(m)
    if len(sort_by_n_connections) == 4:
        central = sort_by_n_connections[3]
        branch1, branch2, branch3 = sorted(
            sort_by_n_connections[:3],
            key=lambda x: top.get_index(x),
        )
        return (
            top.get_index(central),
            top.get_index(branch1),
            top.get_index(branch2),
            top.get_index(branch3),
        )
    return None


def _trim_duplicates(all_matches):
    """Remove redundant sub-graph matches.

    Is there a better way to do this? Like when we format the subgraphs,
    can we impose an ordering so it's easier to eliminate redundant matches?
    """
    trimmed_list = IndexedSet()
    for match in all_matches:
        if (
            match
            and match not in trimmed_list
            and match[::-1] not in trimmed_list
        ):
            trimmed_list.add(match)
    return trimmed_list


def generate_pairs_lists(
    top, molecule=None, sort_key=None, refer_from_scaling_factor=False
):
    """Generate all the pairs lists of the topology or molecular of topology.

    Parameters
    ----------
    top : gmso.Topology
        The Topology where we want to generate the pairs lists from.
    molecule : molecule namedtuple, optional, default=None
        Generate only pairs list of a particular molecule.
    sort_key : function, optional, default=None
        Function used as key for sorting of site pairs. If None is provided
        will used topology.get_index
    refer_from_scaling_factor : bool, optional, default=False
        If True, only generate pair lists of pairs that have a non-zero scaling
        factor value.

    Returns
    -------
    pairs_lists: dict of list
        {"pairs12": pairs12, "pairs13": pairs13, "pairs14": pairs14}

    NOTE: This method assume that the topology has already been loaded with
    angles and dihedrals (through top.identify_connections()). In addition,
    if the refer_from_scaling_factor is True, this method will only generate
    pairs when the corresponding scaling factor is not 0.
    """
    from gmso.external import to_networkx
    from gmso.parameterization.molecule_utils import (
        molecule_angles,
        molecule_bonds,
        molecule_dihedrals,
    )

    nb_scalings, coulombic_scalings = top.scaling_factors

    if sort_key is None:
        sort_key = top.get_index

    graph = to_networkx(top, parse_angles=False, parse_dihedrals=False)

    pairs_dict = dict()
    if refer_from_scaling_factor:
        for i in range(3):
            if nb_scalings[i] or coulombic_scalings[i]:
                pairs_dict[f"pairs1{i+2}"] = list()
    else:
        for i in range(3):
            pairs_dict = {f"pairs1{i+2}": list() for i in range(3)}

    if molecule is None:
        bonds, angles, dihedrals = top.bonds, top.angles, top.dihedrals
    else:
        bonds = molecule_bonds(top, molecule)
        angles = molecule_angles(top, molecule)
        dihedrals = molecule_dihedrals(top, molecule)

    if "pairs12" in pairs_dict:
        for bond in bonds:
            pairs = sorted(bond.connection_members, key=sort_key)
            pairs_dict["pairs12"].append(pairs)

    if "pairs13" in pairs_dict:
        for angle in angles:
            pairs = sorted(
                (angle.connection_members[0], angle.connection_members[-1]),
                key=sort_key,
            )
            if (
                pairs not in pairs_dict["pairs13"]
                and shortest_path_length(graph, pairs[0], pairs[1]) == 2
            ):
                pairs_dict["pairs13"].append(pairs)

    if "pairs14" in pairs_dict:
        for dihedral in dihedrals:
            pairs = sorted(
                (
                    dihedral.connection_members[0],
                    dihedral.connection_members[-1],
                ),
                key=sort_key,
            )
            if (
                pairs not in pairs_dict["pairs14"]
                and shortest_path_length(graph, pairs[0], pairs[1]) == 3
            ):
                pairs_dict["pairs14"].append(pairs)

    for key in pairs_dict:
        pairs_dict[key] = sorted(
            pairs_dict[key],
            key=lambda pairs: (sort_key(pairs[0]), sort_key(pairs[1])),
        )

    return pairs_dict
