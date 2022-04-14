"""TopologyGraph Functions that identify molecules from isomorphism."""
from collections import deque

import networkx as nx


def partition_isomorphic_topology_graphs(graph):
    """Return a dictionary of subgraphs that are partitioned by isomorphism.

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

    subgraphs_gen = (graph.subgraph(c) for c in nx.connected_components(graph))
    subgraphs_list = list(subgraphs_gen)

    graph_queue = deque(subgraphs_list)
    graph_of_interest = graph_queue.popleft()
    isomorphic_elements = {
        graph_of_interest: [],
    }

    last_element_popped = None
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
                graph, graph_of_interest
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
