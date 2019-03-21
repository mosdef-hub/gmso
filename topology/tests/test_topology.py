from copy import deepcopy

import numpy as np
import pytest
import unyt as u
import parmed as pmd

from topology.core.topology import Topology
from topology.core.site import Site
from topology.core.bond import Bond
from topology.core.box import Box
from topology.core.atom_type import AtomType
from topology.external.convert_parmed import from_parmed
from topology.utils.io import get_fn
from topology.testing.utils import allclose
from topology.tests.base_test import BaseTest


class TestTopology(BaseTest):
    def test_new_topology(self):
        top = Topology(name='mytop')
        assert top.name == 'mytop'

    def test_add_site(self):
        top = Topology()
        site = Site(name='site')

        assert top.n_sites == 0
        top.add_site(site)
        assert top.n_sites == 1

    def test_add_connection(self):
        top = Topology()
        site1 = Site(name='site1')
        site2 = Site(name='site2')
        connect = Bond(connection_members=[site1,site2])

        top.add_connection(connect)
        top.add_site(site1)
        top.add_site(site2)


        assert len(top.connection_list) == 1

    def test_add_box(self):
        top = Topology()
        box = Box(2*u.nm*np.ones(3))

        assert top.box is None
        top.box = box
        assert top.box is not None
        assert allclose(top.box.lengths, u.nm*2*np.ones(3))

    def test_positions_dtype(self):
        top = Topology()
        site1 = Site(name='site1')
        top.add_site(site1)

        assert set([type(site.position) for site in top.site_list]) == {u.unyt_array}
        assert set([site.position.units for site in top.site_list]) == {u.nm}

        assert top.positions().dtype == float
        assert top.positions().units == u.nm
        assert isinstance(top.positions(), u.unyt_array)

    def test_eq(self, box, charge):
        top = Topology(name='mytop')

        assert top != box

        diff_name = deepcopy(top)
        diff_name.name = 'othertop'
        assert top != diff_name

        site = Site()
        wrong_n_sites = deepcopy(top)
        assert top == wrong_n_sites
        top.add_site(site)
        assert top != wrong_n_sites

        wrong_position = deepcopy(top)
        top.add_site(Site(position=u.nm*[0, 0, 0]))
        wrong_position.add_site(Site(position=u.nm*[1, 1, 1]))
        assert top != wrong_position

        wrong_charge = deepcopy(top)
        top.add_site(Site(charge=charge))
        wrong_position.add_site(Site(charge=-1*charge))
        assert top != wrong_charge

        wrong_atom_type = deepcopy(top)
        top.add_site(Site(atom_type=AtomType(expression='epsilon*sigma')))
        wrong_atom_type.add_site(Site(atom_type=AtomType(expression='sigma')))
        assert top != wrong_atom_type

        top1 = from_parmed(pmd.load_file(get_fn('ethane.top'),
                                         xyz=get_fn('ethane.gro')))
        top2 = from_parmed(pmd.load_file(get_fn('ethane.top'),
                                         xyz=get_fn('ethane.gro')))
        assert top1 == top2
