from collections import defaultdict


class ConnectionGraph(object):
    """
    A graph that describes the connectivity of sites in a topology

    This code is copied verbatim from mbuild's BondGraph, which is
    designed to mimic the API and partial functionality of NetworkX's
    `Graph` data structure.
     """
    def __init__(self):
        self._adj = defaultdict(set)

    def add_node(self, node):
        if not self.has_node(node):
            self._adj[node] = set()

    def remove_node(self, node):
        adj = self._adj
        for other_node in self.nodes():
            if node in adj[other_node]:
                self.remove_edge(node, other_node)

    def has_node(self, node):
        return node in self._adj

    def nodes(self):
        return [node for node in self._adj]

    def nodes_iter(self):
        for node in self._adj:
            yield node

    def number_of_nodes(self):
        return sum(1 for _ in self.nodes_iter())

    def add_edge(self, node1, node2):
        self._adj[node1].add(node2)
        self._adj[node2].add(node1)

    def remove_edge(self, node1, node2):
        adj = self._adj
        if self.has_node(node1) and self.has_node(node2):
            adj[node1].remove(node2)
            adj[node2].remove(node1)
            if not adj[node1]:
                del adj[node1]
            if not adj[node2]:
                del adj[node2]
        else:
            raise ValueError('There is no edge between {} and {}'.format(
                node1, node2))

    def has_edge(self, node1, node2):
        if self.has_node(node1):
            return node2 in self._adj[node1]

    def edges(self):
        edges = set()
        for node, neighbors in self._adj.items():
            for neighbor in neighbors:
                bond = (node, neighbor) if id(node) < id(neighbor) else (neighbor, node)
                edges.add(bond)
        return list(edges)

    def edges_iter(self):
        for edge in self.edges():
            yield edge

    def number_of_edges(self):
        return sum(1 for _ in self.edges())

    def neighbors(self, node):
        if self.has_node(node):
            return [neighbor for neighbor in self._adj[node]]
        else:
            return []

    def neighbors_iter(self, node):
        if self.has_node(node):
            return (neighbor for neighbor in self._adj[node])
        else:
            return iter(())

    def compose(self, graph):
        adj = self._adj
        for node, neighbors in graph._adj.items():
            if self.has_node(node):
                adj[node].union(neighbors)
            elif neighbors:
                adj[node] = neighbors

    def subgraph(self, nodes):
        new_graph = BondGraph()
        nodes = list(nodes)
        adj = self._adj
        for node in nodes:
            if node not in adj:
                continue
            for neighbor in adj[node]:
                if neighbor in nodes:
                    new_graph.add_edge(node, neighbor)
        return new_graph

def graph_from_top(top):
    """Generate a ConnectionGraph from an existing topology."""
    graph = ConnectionGraph()

    for site in top.sites:
        graph.add_node(site)

    for bond in top.bonds:
        graph.add_edge(
            node1=bond.connection_members[0],
            node2=bond.connection_members[1],
        )

    return graph
