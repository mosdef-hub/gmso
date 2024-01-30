"""Convert to/from NetworkX graphs and GMSO topologies."""

import warnings

import networkx as nx

from gmso.abc.abstract_connection import Connection
from gmso.abc.abstract_site import Site
from gmso.core.bond import Bond
from gmso.core.topology import Topology


def from_networkx(graph):
    """Convert a networkx.Graph to a gmso.Topology.

    Creates a topology from the graph where each node is a site and each
    edge becomes a connection.

    Parameters
    ----------
    graph : networkX.Graph
        networkx.Graph instance that need to be converted

    Returns
    -------
    top : gmso.Topology

    Notes
    -----
    - While a lot of information is lost from converting to a graph object
    (e.g. metadata, mixing rules, etc.), the graph representation is a
    useful way to manipulate and extract connectivity information from
    Topology objects.
    - The edge has a `connection` attribute, which stores the Bond
    object it was created from
    """
    if not isinstance(graph, nx.Graph):
        raise TypeError(
            "Type mismatch, graph object is expected to be "
            "an instance of networkx.Graph, was provided: {}".format(
                type(graph)
            )
        )
    top = Topology()

    node_mapping = dict()

    for node in graph.nodes:
        if not isinstance(node, Site):
            raise TypeError("Nodes must be instances of gmso.abc.Site")
        else:
            top.add_site(node)

    for edge in graph.edges:
        try:
            conn = graph.get_edge_data(*edge)["connection"]
            if isinstance(conn, Connection) and set(edge).issubset(
                set(conn.connection_members)
            ):
                top.add_connection(conn)
        except KeyError:
            conn = Bond(connection_members=edge)
            top.add_connection(conn)

    warnings.simplefilter("once", UserWarning)

    for node in graph.nodes:
        try:
            graph.nodes[node]["angles"] or graph.nodes[node]["dihedrals"]
            warnings.warn("Angle and Dihedral information is not converted.")
        except KeyError:
            pass

    return top


def to_networkx(top, parse_angles=True, parse_dihedrals=True):
    """Convert a gmso.Topology to a networkX.Graph.

    Creates a graph from the topology where each node is a site and each
    edge is a connection.

    Parameters
    ----------
    top : gmso.Topology
        topology.Topology instance that need to be converted
    parse_angles : bool, optional, default=True
        Populate angle field of all nodes
    parse_dihedral : bool, optional default=True
        Populate dihedral field of all nodes

    Returns
    -------
    graph : networkX.Graph

    Notes
    -----
    While a lot of information is lost from converting to a graph object
    (e.g. metadata, mixing rules, etc.), the graph representation is a
    useful way to manipulate and extract connectivity information from
    Topology objects.
    """
    graph = nx.Graph()

    for n in top.sites:
        graph.add_node(n)

    for b in top.bonds:
        graph.add_edge(
            b.connection_members[0], b.connection_members[1], connection=b
        )

    if parse_angles:
        for node in graph.nodes:
            graph.nodes[node]["angles"] = top._get_angles_for(node)

    if parse_dihedrals:
        for node in graph.nodes:
            graph.nodes[node]["dihedrals"] = top._get_dihedrals_for(node)

    return graph
