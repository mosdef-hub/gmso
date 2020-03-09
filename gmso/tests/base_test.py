import pytest
import numpy as np
import mbuild as mb
import unyt as u

from gmso.core.box import Box
from gmso.core.topology import Topology
from gmso.core.element import Hydrogen
from gmso.core.site import Site
from gmso.core.atom_type import AtomType
from gmso.core.forcefield import ForceField
from gmso.external.convert_mbuild import from_mbuild
from gmso.tests.utils import get_path
from gmso.utils.io import get_fn


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

    @pytest.fixture
    def typed_single_ar(self):
        top = from_mbuild(mb.Compound(name="Ar"))

        ff = ForceField(get_fn("ar.xml"))

        for site in top.sites:
            site.atom_type = ff.atom_types["Ar"]

        top.update_topology()
        return top

    @pytest.fixture
    def typed_single_xe_mie(self):
        top = from_mbuild(mb.Compound(name="Xe"))

        ff = ForceField(get_path("noble_mie.xml"))

        for site in top.sites:
            site.atom_type = ff.atom_types["Xe"]

        top.update_topology()
        return top

    @pytest.fixture
    def typed_ar_system(self, ar_system):
        top = ar_system

        ff = ForceField(get_fn('ar.xml'))

        for site in top.sites:
            site.atom_type = ff.atom_types['Ar']

        top.update_topology()

        return top

    @pytest.fixture
    def water_system(self):
        water = mb.load(get_path('tip3p.mol2'))
        water.name = 'water'
        water[0].name = 'opls_111'
        water[1].name = water[2].name = 'opls_112'

        packed_system = mb.fill_box(
                compound=water,
                n_compounds=10,
                box=mb.Box([2, 2, 2])
                )

        return  from_mbuild(packed_system)
