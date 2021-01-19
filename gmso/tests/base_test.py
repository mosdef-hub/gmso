import pytest
import numpy as np
import mbuild as mb
import mbuild.recipes
import unyt as u
import foyer

from gmso.core.box import Box
from gmso.core.topology import Topology
from gmso.core.element import Hydrogen, Oxygen
from gmso.core.atom import Atom
from gmso.core.bond import Bond
from gmso.core.angle import Angle
from gmso.core.atom_type import AtomType
from gmso.core.forcefield import ForceField
from gmso.external import from_mbuild, from_parmed
from gmso.tests.utils import get_path
from gmso.utils.io import get_fn, has_foyer
from gmso.external import from_parmed
from gmso.external.convert_foyer import from_foyer_xml

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
    def ar_system(self, n_ar_system):
        return from_mbuild(n_ar_system())

    @pytest.fixture
    def n_ar_system(self):
        def _topology(n_sites=100):
            ar = mb.Compound(name='Ar')

            packed_system = mb.fill_box(
                compound=ar,
                n_compounds=n_sites,
                box=mb.Box([3, 3, 3]),
            )

            return packed_system

        return _topology

    @pytest.fixture
    def n_typed_xe_mie(self):
        def _typed_topology(n_sites=1):
            xe = mb.Compound(name="Xe")

            packed_system = mb.fill_box(
                compound=xe,
                n_compounds=n_sites,
                box=mb.Box([3, 3, 3]),
            )

            top = from_mbuild(packed_system)

            ff = ForceField(get_path("noble_mie.xml"))

            for site in top.sites:
                site.atom_type = ff.atom_types["Xe"]

            top.update_topology()

            return top

        return _typed_topology

    @pytest.fixture
    def typed_ar_system(self, n_typed_ar_system):
        return n_typed_ar_system()

    @pytest.fixture
    def n_typed_ar_system(self, n_ar_system):
        def _typed_topology(n_sites=100):
            top = from_mbuild(n_ar_system(n_sites=n_sites))

            ff = ForceField(get_fn('ar.xml'))

            for site in top.sites:
                site.atom_type = ff.atom_types['Ar']

            top.update_topology()

            return top

        return _typed_topology

    @pytest.fixture
    def water_system(self):
        water = mb.load(get_path('tip3p.mol2'))
        water.name = 'water'
        water[0].name = 'opls_111'
        water[1].name = water[2].name = 'opls_112'

        packed_system = mb.fill_box(
                compound=water,
                n_compounds=2,
                box=mb.Box([2, 2, 2])
                )

        return  from_mbuild(packed_system)

    @pytest.fixture
    def ethane(self):
        from mbuild.lib.molecules import Ethane
        top = from_mbuild(Ethane())
        return top

    @pytest.fixture
    def typed_ethane(self):
        from mbuild.lib.molecules import Ethane
        mb_ethane = Ethane()
        oplsaa = foyer.Forcefield(name='oplsaa')
        # At this point, we still need to go through
        # parmed Structure, until foyer can perform
        # atomtyping on gmso Topology
        pmd_ethane = oplsaa.apply(mb_ethane)
        top = from_parmed(pmd_ethane)
        top.name = "ethane"
        return top

    @pytest.fixture
    def parmed_ethane(self):
        from mbuild.lib.molecules import Ethane
        compound = Ethane()
        oplsaa = foyer.Forcefield(name='oplsaa')
        pmd_structure = oplsaa.apply(compound)
        return pmd_structure

    @pytest.fixture
    def parmed_methylnitroaniline(self):
        compound = mb.load('CC1=C(C=CC(=C1)[N+](=O)[O-])N', smiles=True)
        oplsaa = foyer.Forcefield(name='oplsaa')
        pmd_structure = oplsaa.apply(compound)
        return pmd_structure

    @pytest.fixture
    def typed_methylnitroaniline(self):
        compound = mb.load('CC1=C(C=CC(=C1)[N+](=O)[O-])N', smiles=True)
        oplsaa = foyer.Forcefield(name='oplsaa')
        pmd_structure = oplsaa.apply(compound)
        top = from_parmed(pmd_structure)
        return top

    @pytest.fixture
    def parmed_chloroethanol(self):
        compound = mb.load('C(CCl)O', smiles=True)
        oplsaa = foyer.Forcefield(name='oplsaa')
        pmd_structure = oplsaa.apply(compound)
        return pmd_structure

    @pytest.fixture
    def typed_chloroethanol(self):
        compound = mb.load('C(CCl)O', smiles=True)
        oplsaa = foyer.Forcefield(name='oplsaa')
        pmd_structure = oplsaa.apply(compound)
        top = from_parmed(pmd_structure)
        return top

    @pytest.fixture
    def parmed_hexane_box(self):
        compound = mb.recipes.Alkane(6)
        compound.name = "HEX"
        compound_box = mb.fill_box(compound, n_compounds=6, box=[6,6,6])
        oplsaa = foyer.Forcefield(name='oplsaa')
        pmd_structure = oplsaa.apply(compound_box, residues="HEX")
        return pmd_structure

    @pytest.fixture
    def typed_water_system(self, water_system):
        top = water_system

        ff = ForceField(get_path('tip3p.xml'))

        element_map = {"O": "opls_111", "H": "opls_112"}

        for atom in top.sites:
            atom.atom_type = ff.atom_types[atom.name]

        for bond in top.bonds:
            bond.bond_type = ff.bond_types["opls_111~opls_112"]

        for subtop in top.subtops:
            angle = Angle(
                connection_members=[site for site in subtop.sites],
                name="opls_112~opls_111~opls_112",
                angle_type=ff.angle_types["opls_112~opls_111~opls_112"]
            )
            top.add_connection(angle)

        top.update_topology()
        return top

    @pytest.fixture
    def foyer_fullerene(self):
       if has_foyer:
          import foyer
          from foyer.tests.utils import get_fn
       from_foyer_xml(get_fn("fullerene.xml"), overwrite=True)
       gmso_ff = ForceField(get_fn("fullerene_gmso.xml"))

       return gmso_ff


    @pytest.fixture
    def foyer_periodic(self):
       if has_foyer:
          import foyer
          from foyer.tests.utils import get_fn
       from_foyer_xml(get_fn("oplsaa-periodic.xml"), overwrite=True)
       gmso_ff = ForceField(get_fn("oplsaa-periodic_gmso.xml"))

       return gmso_ff

    @pytest.fixture
    def foyer_urey_bradley(self):
        if has_foyer:
            import foyer
            from foyer.tests.utils import get_fn
            from_foyer_xml(get_fn("charmm36_cooh.xml"), overwrite=True)
            gmso_ff = ForceField(get_fn("charmm36_cooh_gmso.xml"))

            return gmso_ff

    @pytest.fixture
    def foyer_rb_torsion(self):
        if has_foyer:
            import foyer
            from foyer.tests.utils import get_fn
            from_foyer_xml(get_fn("refs-multi.xml"), overwrite=True, validate_foyer=True)
            gmso_ff = ForceField(get_fn("refs-multi_gmso.xml"))

            return gmso_ff

    @pytest.fixture
    def foyer_periodic_improper(self):
        if has_foyer:
            import foyer
            from foyer.tests.utils import get_fn
            from_foyer_xml(get_fn("improper_dihedral.xml"), overwrite=True, validate_foyer=True)
            gmso_ff = ForceField(get_fn("improper_dihedral_gmso.xml"))

            return gmso_ff

    @pytest.fixture
    def methane(self):
        mytop = Topology()
        c = Atom(name='c')
        h1 = Atom(name='h1')
        h2 = Atom(name='h2')
        h3 = Atom(name='h3')
        h4 = Atom(name='h4')
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
        mytop.update_topology()

        return mytop

    @pytest.fixture
    def ethane(self):
        mytop = Topology()
        c1 = Atom(name='C1')
        h11 = Atom(name='H11')
        h12 = Atom(name='H12')
        h13 = Atom(name='H13')

        c2 = Atom(name='C2')
        h21 = Atom(name='H21')
        h22 = Atom(name='H22')
        h23 = Atom(name='H23')

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
        mytop.update_topology()

        return mytop
