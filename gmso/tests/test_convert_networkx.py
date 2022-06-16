import networkx as nx
import pytest
import unyt as u

import gmso
from gmso.core.atom import Atom
from gmso.core.topology import Topology as Top
from gmso.external.convert_networkx import from_networkx, to_networkx
from gmso.tests.base_test import BaseTest


class TestConvertNetworkX(BaseTest):
    def test_to_networkx_ethane(self, ethane_from_scratch):
        ethane_to_nx = to_networkx(ethane_from_scratch)

        assert ethane_from_scratch.n_sites == ethane_to_nx.number_of_nodes()
        assert ethane_from_scratch.n_bonds == ethane_to_nx.number_of_edges()

        assert set(ethane_from_scratch.sites) == set(ethane_to_nx.nodes)
        for site in ethane_from_scratch.sites:
            assert set(ethane_to_nx.nodes[site]["angles"]) == set(
                ethane_from_scratch._get_angles_for(site)
            )
            assert set(ethane_to_nx.nodes[site]["dihedrals"]) == set(
                ethane_from_scratch._get_dihedrals_for(site)
            )

    def test_from_networkx_ethane(self, ethane_from_scratch):
        ethane_to_nx = to_networkx(ethane_from_scratch)
        ethane_from_nx = from_networkx(ethane_to_nx)

        assert ethane_from_scratch.n_sites == ethane_from_nx.n_sites
        assert ethane_from_scratch.n_bonds == ethane_from_nx.n_bonds

        assert set(ethane_from_scratch.sites) == set(ethane_from_nx.sites)
        assert set(ethane_from_scratch.bonds) == set(ethane_from_nx.bonds)

    def test_from_networkx_argon(self, ar_system):
        ar_to_nx = to_networkx(ar_system)
        ar_from_nx = from_networkx(ar_to_nx)

        assert ar_system.n_sites == ar_from_nx.n_sites
        assert ar_system.n_bonds == ar_from_nx.n_bonds

        assert set(ar_system.sites) == set(ar_from_nx.sites)
        assert set(ar_system.bonds) == set(ar_from_nx.bonds)

    def test_from_networkx_not_graph_object(self):
        fake_graph = [23, 23]
        with pytest.raises(TypeError):
            from_networkx(fake_graph)

    def test_from_networkx_water_box(self, water_system):
        water_to_nx = to_networkx(water_system)
        water_from_nx = from_networkx(water_to_nx)

        assert water_system.n_sites == water_from_nx.n_sites
        assert water_system.n_bonds == water_from_nx.n_bonds

        assert set(water_system.sites) == set(water_from_nx.sites)
        assert set(water_system.bonds) == set(water_from_nx.bonds)

        assert nx.number_connected_components(water_to_nx) == len(
            water_system.unique_site_labels(
                label_type="molecule", name_only=False
            )
        )

    def test_from_networkx_without_connections(self):
        g = nx.Graph()
        g.add_edge(Atom(), Atom())
        top = from_networkx(g)
        assert top.n_connections == 1

    def test_from_networkx_arbitrary_graph(self):
        test_graph = nx.grid_2d_graph(2, 2)
        with pytest.raises(TypeError) as e:
            from_networkx(test_graph)
