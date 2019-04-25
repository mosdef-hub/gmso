from topology.core.topology import Topology
from topology.core.site import Site
from topology.core.subtopology import SubTopology

from topology.tests.base_test import BaseTest


class TestBond(BaseTest):
    def test_subtop_add_site(self):
        subtop = SubTopology()
        site = Site()

        subtop.add_site(site)

        assert subtop.n_sites == 1

    def test_subtop_assign_parent(self):
        subtop = SubTopology()
        top = Topology()

        subtop.parent = top

        assert subtop.parent is not None
