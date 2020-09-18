import numpy as np
import unyt as u

from gmso.core.topology import Topology
from gmso.core.atom import Atom
from gmso.core.subtopology import SubTopology
from gmso.tests.base_test import BaseTest


class TestSubTopology(BaseTest):
    def test_subtop_init(self):
        default = SubTopology()
        assert default.name == 'Sub-Topology'
        assert default.parent is None
        assert default.n_sites == 0

        named = SubTopology(name='CoolSubTopology')
        assert named.name == 'CoolSubTopology'
        assert named.parent is None
        assert named.n_sites == 0

        myparent = Topology()
        with_parent = SubTopology(parent=myparent)
        assert with_parent.name == 'Sub-Topology'
        assert with_parent.parent == myparent
        assert with_parent.n_sites == 0

    def test_subtop_setters(self):
        subtop = SubTopology()
        assert subtop.name == 'Sub-Topology'
        assert subtop.parent is None
        assert subtop.n_sites == 0

        subtop.name = 'NewSubTopology'
        assert subtop.name == 'NewSubTopology'

        newparent = Topology()
        subtop.parent = newparent
        assert subtop.parent == newparent

    def test_subtop_add_site(self):
        subtop = SubTopology()
        site = Atom()

        subtop.add_site(site)

        assert subtop.n_sites == 1

    def test_subtop_assign_parent(self):
        subtop = SubTopology()
        top = Topology()

        subtop.parent = top

        assert subtop.parent is not None

    def test_subtop_add_site_parent(self):
        top = Topology()
        subtop = SubTopology(parent=top)
        site = Atom()

        subtop.add_site(site)

        assert subtop.n_sites == 1
        assert subtop.parent.n_sites == 1
        assert top.n_sites == 1

    def test_add_site_parent(self):
        top = Topology()
        subtop = SubTopology()
        site1 = Atom(position=u.nm * np.zeros(3))
        site2 = Atom(position=u.nm * np.ones(3))
        top.add_subtopology(subtop)

        assert top.n_sites == 0
        assert subtop.n_sites == 0
        subtop.add_site(site1)
        assert top.n_sites == 1
        assert subtop.n_sites == 1
        top.add_site(site2)
        assert top.n_sites == 2
        assert subtop.n_sites == 1
