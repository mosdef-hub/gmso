import networkx as nx

from gmso.core.bond import Bond
from gmso.core.angle import Angle


def identify_connections(top):
    """ Identify all possible connections within a topology

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

    angle_matches = _detect_angles(compound_line_graph)
    # Todo: Uncomment when dihedral and improper class get implemented
    #dihedral_matches = _detect_dihedrals(compound_line_graph)
    #improper_matches = _detect_impropers(compound_line_graph)

    _add_angles(top, angle_matches)
    #_add_dihedral(top, dihedral_matches)
    #_add_impropers(top, improper_matches)

    return top

def _detect_angles(compound_line_graph):
    angle = nx.Graph()
    angle.add_edge(0, 1)

    matcher = nx.algorithms.isomorphism.GraphMatcher(compound_line_graph, angle)

    angle_matches = []
    for m in matcher.subgraph_isomorphisms_iter():
        new_connection = _format_subgraph_angle(m)
        angle_matches.append(new_connection)
    angle_matches = _trim_duplicates(angle_matches)

    return angle_matches

def _add_angles(top, angle_matches):
    for angle_tuple in angle_matches:
        to_add_angle = Angle(connection_members=[*angle_tuple])
        top.add_connection(to_add_angle, update_types=False)

def _add_dihedrals(top, dihedral_matches):
    for dihedral_tuple in dihedral_matches:
        to_add_dihedral = Dihedral(connection_members=[*dihedral_tuple])
        top.add_connection(to_add_dihedral, update_types=False)

def _add_impropers(top, improper_matches):
    for improper_tuple in improper_matches:
        to_add_improper = Improper(connection_members=[*improper_tuple])
        top.add_connection(to_add_improper, update_types=False)

def _detect_dihedrals(compound_line_graph):
    dihedral = nx.Graph()
    dihedral.add_edge(0,1)
    dihedral.add_edge(1,2)

    matcher = nx.algorithms.isomorphism.GraphMatcher(compound_line_graph, dihedral)

    dihedral_matches = []
    for m in matcher.subgraph_isomorphisms_iter():
        new_connection = _format_subgraph_dihedral(m)
        dihedral_matches.append(new_connection)
    dihedral_matches = _trim_duplicates(dihedral_matches)

    return dihedral_matches

def _detect_impropers(compound_line_graph):
    improper = nx.Graph()
    improper.add_edge(0,1)
    improper.add_edge(1,2)
    improper.add_edge(0,2)

    matcher = nx.algorithms.isomorphism.GraphMatcher(compound_line_graph, improper)

    improper_matches = []
    for m in matcher.subgraph_isomorphisms_iter():
        new_connection = _format_subgraph_improper(m)
        if new_connection is not None:
            improper_matches.append(new_connection)
    improper_matches = _trim_duplicates(improper_matches)

    return improper_matches

def _format_subgraph_angle(m):
    """ Format the angle subgraph

    Since we are matching compound line graphs,
    back out the actual nodes, not just the edges

    Parameters
    ----------
    m : dict
        keys are the compound line graph nodes
        Values are the sub-graph matches (to the angle, dihedral, or improper)

    Returns
    ------
    connection : list of nodes, in order of bonding
        (start, middle, end) """

    small = nx.Graph()
    for k, v in m.items():
        small.add_edge(k[0],k[1])
    sort_by_n_connections = sorted(small.adj, key=lambda x:len(small[x]))
    start = sort_by_n_connections[0]
    end = sort_by_n_connections[1]
    middle = sort_by_n_connections[2]

    return [start, middle, end]

def _format_subgraph_dihedral(m):
    """ Format the dihedral subgraph

    Since we are matching compound line graphs,
    back out the actual nodes, not just the edges

    Parameters
    ----------
    m : dict
        keys are the compound line graph nodes
        Values are the sub-graph matches (to the angle, dihedral, or improper)

    Returns
    ------
    connection : list of nodes, in order of bonding
        (start, mid1, mid2, end)
        """
    small = nx.Graph()
    for k,v in m.items():
        small.add_edge(k[0], k[1])
    sort_by_n_connections = sorted(small.adj, key=lambda x:len(small[x]))
    start = sort_by_n_connections[0]
    if sort_by_n_connections[2] in small.neighbors(start):
        mid1 = sort_by_n_connections[2]
        mid2 = sort_by_n_connections[3]
    else:
        mid1 = sort_by_n_connections[3]
        mid2 = sort_by_n_connections[2]

    end = sort_by_n_connections[1]
    return [start,mid1, mid2, end]

def _format_subgraph_improper(m):
    """ Format the dihedral subgraph

    Since we are matching compound line graphs,
    back out the actual nodes, not just the edges

    Parameters
    ----------
    m : dict
        keys are the compound line graph nodes
        Values are the sub-graph matches (to the angle, dihedral, or improper)

    Returns
    ------
    connection : list of nodes, in order of bonding
        (central, branch1, branch2, branch3)

    Notes
    ------
    Given the way impropers are matched, sometimes a cyclic 3-ring system gets returned
    """
    small = nx.Graph()
    for k,v in m.items():
        small.add_edge(k[0], k[1])
    sort_by_n_connections = sorted(small.adj, key=lambda x:len(small[x]))
    if len(sort_by_n_connections) == 4:
        central = sort_by_n_connections[3]
        branch1, branch2, branch3 = sorted(sort_by_n_connections[:3])
        return [central, branch1, branch2, branch3]
    return None

def _trim_duplicates(all_matches):
    """ Remove redundant sub-graph matches

    Is there a better way to do this? Like when we format the subgraphs,
    can we impose an ordering so it's easier to eliminate redundant matches?
    """
    trimmed_list = []
    for match in all_matches:
        if match not in trimmed_list and match[::-1] not in trimmed_list:
            trimmed_list.append(match)
    return trimmed_list

