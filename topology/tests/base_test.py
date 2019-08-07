import pytest
import numpy as np
import unyt as u

from topology.core.box import Box
from topology.core.topology import Topology
from topology.core.element import Element
from topology.core.site import Site
from topology.core.atom_type import AtomType
from topology.core.bond import Bond


class BaseTest:
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def lengths(self):
        return u.nm * np.ones(3)

    @pytest.fixture
    def angles(self):
        return u.degree * [90, 90, 90]

    @pytest.fixture
    def charge(self):
        return u.elementary_charge * 1

    @pytest.fixture
    def mass(self):
        return 1 * u.gram/u.mol

    @pytest.fixture
    def box(self):
        return Box(lengths=u.nm*np.ones(3))

    @pytest.fixture
    def top(self):
        return Topology(name='mytop')

    @pytest.fixture
    def topology_site(self):
        def _topology(sites=1):
            top = Topology()
            top.box = Box(lengths=[1, 1, 1])
            H = Element(name='H', symbol='H', mass=1)
            site1 = Site(name='site1',
                         element=H,
                         atom_type=AtomType(name="at1",
                                            mass=H.mass)
                         )
            for i in range(sites):
                top.add_site(site1)

            return top

        return _topology

    @pytest.fixture
    def methane(self):
        mytop = Topology()
        c = Site(name='c')
        h1 = Site(name='h1')
        h2 = Site(name='h2')
        h3 = Site(name='h3')
        h4 = Site(name='h4')
        ch1 = Bond(connection_members=[c,h1])
        ch2 = Bond(connection_members=[c,h2])
        ch3 = Bond(connection_members=[c,h3])
        ch4 = Bond(connection_members=[c,h4])
        mytop.add_site(c, update_types=False)
        mytop.add_site(h1, update_types=False)
        mytop.add_site(h2, update_types=False)
        mytop.add_site(h3, update_types=False)
        mytop.add_site(h4, update_types=False)
        mytop.add_connection(ch1, update_types=False)
        mytop.add_connection(ch2, update_types=False)
        mytop.add_connection(ch3, update_types=False)
        mytop.add_connection(ch4, update_types=False)
        mytop.update_top(update_types=False)

        return mytop

    @pytest.fixture
    def ethane(self):
        mytop = Topology()
        c1 = Site(name='C1')
        h11 = Site(name='H11')
        h12 = Site(name='H12')
        h13 = Site(name='H13')

        c2 = Site(name='C2')
        h21 = Site(name='H21')
        h22 = Site(name='H22')
        h23 = Site(name='H23')

        c1h11 = Bond(connection_members=[c1, h11])
        c1h12 = Bond(connection_members=[c1, h12])
        c1h13 = Bond(connection_members=[c1, h13])

        c2h21 = Bond(connection_members=[c2, h21])
        c2h22 = Bond(connection_members=[c2, h22])
        c2h23 = Bond(connection_members=[c2, h23])

        c1c2 = Bond(connection_members=[c1,c2])

        mytop.add_connection(c1h11, update_types=False)
        mytop.add_connection(c1h12, update_types=False)
        mytop.add_connection(c1h13, update_types=False)

        mytop.add_connection(c2h21, update_types=False)
        mytop.add_connection(c2h22, update_types=False)
        mytop.add_connection(c2h23, update_types=False)

        mytop.add_connection(c1c2, update_types=False)
        mytop.update_top(update_types=False)

        return mytop
