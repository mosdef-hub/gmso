"""Convert to/from NetworkX graphs and GMSO topologies."""

import logging

import networkx as nx

from gmso.abc.abstract_connection import Connection
from gmso.abc.abstract_site import Site
from gmso.core.bond import Bond
from gmso.core.topology import Topology

logger = logging.getLogger(__name__)


def from_networkx(graph: nx.Graph) -> Topology:
    """Convert a ``networkx.Graph`` to a :class:`~gmso.Topology`.

    Each graph node must be a :class:`~gmso.abc.abstract_site.Site` and
    each edge becomes a :class:`~gmso.Bond`.

    Parameters
    ----------
    graph : networkx.Graph
        Graph whose nodes are :class:`~gmso.abc.abstract_site.Site` instances.
        Edge data may optionally contain a ``'connection'`` key holding a
        :class:`~gmso.abc.abstract_connection.Connection` object.

    Returns
    -------
    gmso.Topology
        Topology containing the sites and bonds from *graph*.

    Notes
    -----
    Much information is lost when converting a topology to a graph (mixing
    rules, metadata, etc.).  The edge ``'connection'`` attribute stores the
    original :class:`~gmso.Bond` and is used to reconstruct it during
    the reverse conversion.

    Angle and dihedral information stored as node attributes is detected but
    not converted; a warning is logged if such data is found.

    Raises
    ------
    TypeError
        If *graph* is not a ``networkx.Graph`` or if any node is not a
        :class:`~gmso.abc.abstract_site.Site`.
    """
    if not isinstance(graph, nx.Graph):
        raise TypeError(
            "Type mismatch, graph object is expected to be "
            "an instance of networkx.Graph, was provided: {}".format(type(graph))
        )
    top = Topology()

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

    for node in graph.nodes:
        try:
            graph.nodes[node]["angles"] or graph.nodes[node]["dihedrals"]
            logger.info("Angle and Dihedral information is not converted.")
        except KeyError:
            pass

    return top


def to_networkx(
    top: Topology,
    parse_angles: bool = True,
    parse_dihedrals: bool = True,
) -> nx.Graph:
    """Convert a :class:`~gmso.Topology` to a ``networkx.Graph``.

    Each site becomes a graph node and each bond becomes an edge.

    Parameters
    ----------
    top : gmso.Topology
        The topology to convert.
    parse_angles : bool, optional, default=True
        When ``True``, populate an ``'angles'`` attribute on every node
        with the list of angles that include that site.
    parse_dihedrals : bool, optional, default=True
        When ``True``, populate a ``'dihedrals'`` attribute on every node
        with the list of dihedrals that include that site.

    Returns
    -------
    networkx.Graph
        Graph where nodes are :class:`~gmso.Atom` objects and edges carry
        a ``'connection'`` attribute holding the corresponding
        :class:`~gmso.Bond`.

    Notes
    -----
    Converting to a graph loses metadata, mixing rules, and all potential
    type information.  The graph is primarily useful for connectivity
    analysis.
    """
    graph = nx.Graph()

    for n in top.sites:
        graph.add_node(n)

    for b in top.bonds:
        graph.add_edge(b.connection_members[0], b.connection_members[1], connection=b)

    if parse_angles:
        for node in graph.nodes:
            graph.nodes[node]["angles"] = top._get_angles_for(node)

    if parse_dihedrals:
        for node in graph.nodes:
            graph.nodes[node]["dihedrals"] = top._get_dihedrals_for(node)

    return graph
