"""TopologyGraph Functions that identify molecules from isomorphism."""

from collections import deque

import networkx as nx


def top_node_match(n1, n2):
    """Match two nodes in the topology graph based on their elements."""
    return n1["atom_data"].element == n2["atom_data"].element


def partition_isomorphic_topology_graphs(graph):
    """Return a collection of isomorphic sets of the subgraphs of the Topology Graph.

    Parameters
    ----------
    graph: foyer.topology_graph.TopologyGraph
        The networkx subclassed TopologyGraph with data identifying the nodes
        and edges that make up a topology atom and bonds structure

    Returns
    -------
    isomorphic_elements: dict
        The keys are unique disconnected graphs, and the values are all
        identical subgraphs in the graph

    Notes
    -----
    See https://github.com/networkx/networkx/blob/main/networkx/algorithms/isomorphism/isomorphvf2.py
    from the networkx documentation about identifying isomorphic components
    """
    graph_queue = deque(
        graph.subgraph(c) for c in nx.connected_components(graph)
    )

    graph_of_interest = graph_queue.popleft()
    isomorphic_elements = {
        graph_of_interest: [],
    }

    count = 0
    first_mismatch = None

    while graph_queue:
        if graph_queue[0] == first_mismatch:
            count = 0
            graph_of_interest = graph_queue.popleft()
            isomorphic_elements[graph_of_interest] = []
        if graph_queue:
            graph = graph_queue.popleft()
            matcher = nx.algorithms.isomorphism.GraphMatcher(
                graph, graph_of_interest, node_match=top_node_match
            )
            if matcher.is_isomorphic():
                isomorphic_elements[graph_of_interest].append(
                    (graph, matcher.mapping)
                )
            else:
                if count == 0:
                    first_mismatch = graph
                graph_queue.append(graph)
                count += 1
    return isomorphic_elements
