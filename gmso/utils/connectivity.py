"""Module supporting various connectivity methods and operations."""
import networkx as nx

from gmso.core.angle import Angle
from gmso.core.dihedral import Dihedral
from gmso.core.improper import Improper

CONNS = {"angle": Angle, "dihedral": Dihedral, "improper": Improper}

EDGES = {
    "angle": ((0, 1),),
    "dihedral": ((0, 1), (1, 2)),
    "improper": ((0, 1), (0, 2), (1, 2)),
}


def identify_connections(top):
    """Identify all possible connections within a topology.

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

    angle_matches = _detect_connections(compound_line_graph, type_="angle")
    dihedral_matches = _detect_connections(
        compound_line_graph, type_="dihedral"
    )
    improper_matches = _detect_connections(
        compound_line_graph, type_="improper"
    )

    for conn_matches, conn_type in zip(
        (angle_matches, dihedral_matches, improper_matches),
        ("angle", "dihedral", "improper"),
    ):
        if conn_matches:
            _add_connections(top, conn_matches, conn_type=conn_type)

    return top


def _add_connections(top, matches, conn_type):
    """Add connections to the topology."""
    for tuple_ in matches:
        to_add_conn = CONNS[conn_type](connection_members=[*tuple_])
        top.add_connection(to_add_conn, update_types=False)


def _detect_connections(compound_line_graph, type_="angle"):
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

    conn_matches = []
    for m in matcher.subgraph_isomorphisms_iter():
        new_connection = formatter_fns[type_](m)
        conn_matches.append(new_connection)

    if conn_matches:
        conn_matches = _trim_duplicates(conn_matches)

    return conn_matches


def _get_sorted_by_n_connections(m):
    """Return sorted by n connections for the matching graph."""
    small = nx.Graph()
    for k, v in m.items():
        small.add_edge(k[0], k[1])
    return sorted(small.adj, key=lambda x: len(small[x])), small


def _format_subgraph_angle(m):
    """Format the angle subgraph.

    Since we are matching compound line graphs,
    back out the actual nodes, not just the edges

    Parameters
    ----------
    m : dict
        keys are the compound line graph nodes
        Values are the sub-graph matches (to the angle, dihedral, or improper)

    Returns
    -------
    connection : list of nodes, in order of bonding
        (start, middle, end)
    """
    (sort_by_n_connections, _) = _get_sorted_by_n_connections(m)
    start = sort_by_n_connections[0]
    end = sort_by_n_connections[1]
    middle = sort_by_n_connections[2]
    return [start, middle, end]


def _format_subgraph_dihedral(m):
    """Format the dihedral subgraph.

    Since we are matching compound line graphs,
    back out the actual nodes, not just the edges

    Parameters
    ----------
    m : dict
        keys are the compound line graph nodes
        Values are the sub-graph matches (to the angle, dihedral, or improper)

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
    return [start, mid1, mid2, end]


def _format_subgraph_improper(m):
    """Format the improper dihedral subgraph.

    Since we are matching compound line graphs,
    back out the actual nodes, not just the edges

    Parameters
    ----------
    m : dict
        keys are the compound line graph nodes
        Values are the sub-graph matches (to the angle, dihedral, or improper)

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
        branch1, branch2, branch3 = sorted(sort_by_n_connections[:3])
        return [central, branch1, branch2, branch3]
    return None


def _trim_duplicates(all_matches):
    """Remove redundant sub-graph matches.

    Is there a better way to do this? Like when we format the subgraphs,
    can we impose an ordering so it's easier to eliminate redundant matches?
    """
    trimmed_list = []
    for match in all_matches:
        if (
            match
            and match not in trimmed_list
            and match[::-1] not in trimmed_list
        ):
            trimmed_list.append(match)
    return trimmed_list
