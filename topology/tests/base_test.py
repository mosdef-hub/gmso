import pytest
import numpy as np
import mbuild as mb
import unyt as u

from topology.core.box import Box
from topology.core.topology import Topology
from topology.core.element import Hydrogen
from topology.core.site import Site
from topology.core.atom_type import AtomType
from topology.external.convert_mbuild import from_mbuild


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
            H = Hydrogen
            site1 = Site(name='site1',
                         element=H,
                         atom_type=AtomType(name="at1",
                                            mass=H.mass),
                         )
            for i in range(sites):
                top.add_site(site1)

            return top

        return _topology


    @pytest.fixture
    def ar_system(self):
        ar = mb.Compound(name='Ar')

        packed_system = mb.fill_box(
            compound=ar,
            n_compounds=100,
            box=mb.Box([3, 3, 3]),
        )

        return from_mbuild(packed_system)
