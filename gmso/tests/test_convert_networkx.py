import pytest
import unyt as u

import gmso
from gmso.core.topology import Topology as Top
from gmso.core.subtopology import SubTopology as SubTop
from gmso.core.atom import Atom
from gmso.external.convert_networkx import from_networkx, to_networkx
from gmso.tests.base_test import BaseTest

class TestConvertNetworkX(BaseTest):
    def test_to_networkx_ethane(self, ethane):
        ethane_to_nx = to_networkx(ethane)

        assert ethane.n_sites == ethane_to_nx.number_of_nodes()
        assert ethane.n_bonds == ethane_to_nx.number_of_edges()

        assert set(ethane.sites) == set(ethane_to_nx.nodes)

    def test_from_networkx_ethane(self, ethane):
        ethane_to_nx = to_networkx(ethane)
        ethane_from_nx = from_networkx(ethane_to_nx)

        assert ethane.n_sites == ethane_from_nx.n_sites
        assert ethane.n_bonds == ethane_from_nx.n_bonds

        assert set(ethane.sites) == set(ethane_from_nx.sites)
        assert set(ethane.bonds) == set(ethane_from_nx.bonds)

