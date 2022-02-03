import foyer
import mbuild as mb
import mbuild.recipes
import numpy as np
import pytest
import unyt as u

from gmso.core.angle import Angle
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.box import Box
from gmso.core.dihedral import Dihedral
from gmso.core.element import Hydrogen, Oxygen
from gmso.core.forcefield import ForceField
from gmso.core.improper import Improper
from gmso.core.pairpotential_type import PairPotentialType
from gmso.core.topology import Topology
from gmso.external import from_mbuild, from_parmed
from gmso.external.convert_foyer_xml import from_foyer_xml
from gmso.tests.utils import get_path
from gmso.utils.io import get_fn, has_foyer


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
        return 1 * u.gram / u.mol

    @pytest.fixture
    def box(self):
        return Box(lengths=u.nm * np.ones(3))

    @pytest.fixture
    def top(self):
        return Topology(name="mytop")

    @pytest.fixture
    def ar_system(self, n_ar_system):
        return from_mbuild(n_ar_system())

    @pytest.fixture
    def n_ar_system(self):
        def _topology(n_sites=100):
            ar = mb.Compound(name="Ar")

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

            ff = ForceField(get_fn("ar.xml"))

            for site in top.sites:
                site.atom_type = ff.atom_types["Ar"]

            top.update_topology()

            return top

        return _typed_topology

    @pytest.fixture
    def water_system(self):
        water = mb.load(get_path("tip3p.mol2"))
        water.name = "water"
        water[0].name = "opls_111"
        water[1].name = water[2].name = "opls_112"

        packed_system = mb.fill_box(
            compound=water, n_compounds=2, box=mb.Box([2, 2, 2])
        )

        return from_mbuild(packed_system)

    @pytest.fixture
    def ethane(self):
        from mbuild.lib.molecules import Ethane

        top = from_mbuild(Ethane())
        return top

    @pytest.fixture
    def typed_ethane(self):
        from mbuild.lib.molecules import Ethane

        mb_ethane = Ethane()
        oplsaa = foyer.Forcefield(name="oplsaa")
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
        oplsaa = foyer.Forcefield(name="oplsaa")
        pmd_structure = oplsaa.apply(compound)
        return pmd_structure

    @pytest.fixture
    def parmed_methylnitroaniline(self):
        compound = mb.load("CC1=C(C=CC(=C1)[N+](=O)[O-])N", smiles=True)
        oplsaa = foyer.Forcefield(name="oplsaa")
        pmd_structure = oplsaa.apply(compound)
        return pmd_structure

    @pytest.fixture
    def typed_methylnitroaniline(self):
        compound = mb.load("CC1=C(C=CC(=C1)[N+](=O)[O-])N", smiles=True)
        oplsaa = foyer.Forcefield(name="oplsaa")
        pmd_structure = oplsaa.apply(compound)
        top = from_parmed(pmd_structure)
        return top

    @pytest.fixture
    def parmed_chloroethanol(self):
        compound = mb.load("C(CCl)O", smiles=True)
        oplsaa = foyer.Forcefield(name="oplsaa")
        pmd_structure = oplsaa.apply(compound)
        return pmd_structure

    @pytest.fixture
    def typed_chloroethanol(self):
        compound = mb.load("C(CCl)O", smiles=True)
        oplsaa = foyer.Forcefield(name="oplsaa")
        pmd_structure = oplsaa.apply(compound)
        top = from_parmed(pmd_structure)
        return top

    @pytest.fixture
    def parmed_hexane_box(self):
        compound = mb.recipes.Alkane(6)
        compound.name = "HEX"
        compound_box = mb.fill_box(compound, n_compounds=6, box=[6, 6, 6])
        oplsaa = foyer.Forcefield(name="oplsaa")
        pmd_structure = oplsaa.apply(compound_box, residues="HEX")
        return pmd_structure

    @pytest.fixture
    def typed_water_system(self, water_system):
        top = water_system

        ff = ForceField(get_path("tip3p.xml"))

        element_map = {"O": "opls_111", "H": "opls_112"}

        for atom in top.sites:
            atom.atom_type = ff.atom_types[atom.name]

        for bond in top.bonds:
            bond.bond_type = ff.bond_types["opls_111~opls_112"]

        for subtop in top.subtops:
            angle = Angle(
                connection_members=[site for site in subtop.sites],
                name="opls_112~opls_111~opls_112",
                angle_type=ff.angle_types["opls_112~opls_111~opls_112"],
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
        gmso_ff = ForceField("fullerene_gmso.xml")

        return gmso_ff

    @pytest.fixture
    def foyer_periodic(self):
        if has_foyer:
            import foyer
            from foyer.tests.utils import get_fn
        from_foyer_xml(get_fn("oplsaa-periodic.xml"), overwrite=True)
        gmso_ff = ForceField("oplsaa-periodic_gmso.xml")

        return gmso_ff

    @pytest.fixture
    def foyer_urey_bradley(self):
        if has_foyer:
            import foyer
            from foyer.tests.utils import get_fn

            from_foyer_xml(get_fn("charmm36_cooh.xml"), overwrite=True)
            gmso_ff = ForceField("charmm36_cooh_gmso.xml")

            return gmso_ff

    @pytest.fixture
    def foyer_rb_torsion(self):
        if has_foyer:
            import foyer
            from foyer.tests.utils import get_fn

            from_foyer_xml(
                get_fn("refs-multi.xml"), overwrite=True, validate_foyer=True
            )
            gmso_ff = ForceField("refs-multi_gmso.xml")

            return gmso_ff

    @pytest.fixture
    def methane(self):
        mytop = Topology()
        c = Atom(name="c")
        h1 = Atom(name="h1")
        h2 = Atom(name="h2")
        h3 = Atom(name="h3")
        h4 = Atom(name="h4")
        ch1 = Bond(connection_members=[c, h1])
        ch2 = Bond(connection_members=[c, h2])
        ch3 = Bond(connection_members=[c, h3])
        ch4 = Bond(connection_members=[c, h4])
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
    def ethane_from_scratch(self):
        mytop = Topology()
        c1 = Atom(name="C1")
        h11 = Atom(name="H11")
        h12 = Atom(name="H12")
        h13 = Atom(name="H13")

        c2 = Atom(name="C2")
        h21 = Atom(name="H21")
        h22 = Atom(name="H22")
        h23 = Atom(name="H23")

        c1h11 = Bond(connection_members=[c1, h11])
        c1h12 = Bond(connection_members=[c1, h12])
        c1h13 = Bond(connection_members=[c1, h13])

        c2h21 = Bond(connection_members=[c2, h21])
        c2h22 = Bond(connection_members=[c2, h22])
        c2h23 = Bond(connection_members=[c2, h23])

        c1c2 = Bond(connection_members=[c1, c2])

        mytop.add_connection(c1h11, update_types=False)
        mytop.add_connection(c1h12, update_types=False)
        mytop.add_connection(c1h13, update_types=False)

        mytop.add_connection(c2h21, update_types=False)
        mytop.add_connection(c2h22, update_types=False)
        mytop.add_connection(c2h23, update_types=False)

        mytop.add_connection(c1c2, update_types=False)
        mytop.update_topology()

        return mytop

    @pytest.fixture(scope="session")
    def are_equivalent_atoms(self):
        def test_atom_equality(atom1, atom2):
            if not all(isinstance(x, Atom) for x in [atom1, atom2]):
                return False
            if atom1 is atom2:
                return True
            equal = (
                lambda x1, x2: u.allclose_units(x1, x2)
                if isinstance(x1, u.unyt_array) and isinstance(x2, u.unyt_array)
                else x1 == x2
            )
            for prop in atom1.dict(by_alias=True):
                if not equal(atom2.dict().get(prop), atom1.dict().get(prop)):
                    return False
            return True

        return test_atom_equality

    @pytest.fixture(scope="session")
    def are_equivalent_connections(self, are_equivalent_atoms):
        connection_types_attrs_map = {
            Bond: "bond_type",
            Angle: "angle_type",
            Dihedral: "dihedral_type",
            Improper: "improper_type",
        }

        def test_connection_equality(conn1, conn2):
            if not type(conn1) == type(conn2):
                return False
            conn1_eq_members = conn1.equivalent_members()
            conn2_eq_members = conn2.equivalent_members()
            have_eq_members = False
            for conn2_eq_member_tuple in conn2_eq_members:
                for conn1_eq_member_tuple in conn1_eq_members:
                    if any(
                        are_equivalent_atoms(member1, member2)
                        for member1, member2 in zip(
                            conn1_eq_member_tuple, conn2_eq_member_tuple
                        )
                    ):
                        have_eq_members = True

            if not have_eq_members:
                return False
            if conn1.name != conn2.name:
                return False
            if getattr(
                conn1, connection_types_attrs_map[type(conn1)]
            ) != getattr(conn2, connection_types_attrs_map[type(conn2)]):
                return False
            return True

        return test_connection_equality

    @pytest.fixture(scope="session")
    def have_equivalent_boxes(self):
        def test_box_equivalence(top1, top2):
            if top1.box and top2.box:
                return u.allclose_units(
                    top1.box.lengths, top2.box.lengths
                ) and u.allclose_units(top1.box.angles, top2.box.angles)
            elif not top1.box and not top2.box:
                return True
            else:
                return False

        return test_box_equivalence

    @pytest.fixture(scope="session")
    def are_equivalent_topologies(
        self,
        have_equivalent_boxes,
        are_equivalent_atoms,
        are_equivalent_connections,
    ):
        def test_topology_equivalence(top1, top2):
            if top1.n_sites != top2.n_sites:
                return False, "Unequal number of sites"
            if top1.n_bonds != top2.n_bonds:
                return False, "Unequal number of bonds"
            if top1.n_angles != top2.n_angles:
                return False, "Unequal number of angles"
            if top1.n_dihedrals != top2.n_dihedrals:
                return False, "Unequal number of dihedrals"
            if top1.n_impropers != top2.n_impropers:
                return False, "Unequal number of impropers"
            if top1.name != top2.name:
                return False, "Dissimilar names"

            if top1.scaling_factors != top2.scaling_factors:
                return False, f"Mismatch in scaling factors"

            if not have_equivalent_boxes(top1, top2):
                return (
                    False,
                    "Non equivalent boxes, differing in lengths and angles",
                )

            for atom1, atom2 in zip(top1.sites, top2.sites):
                if not are_equivalent_atoms(atom1, atom2):
                    return False, f"Non equivalent atoms {atom1, atom2}"

            # Note: In these zipped iterators, index matches are implicitly checked
            for bond1, bond2 in zip(top1.bonds, top2.bonds):
                if not are_equivalent_connections(bond1, bond2):
                    return False, f"Non equivalent bonds {bond1, bond2}"

            for angle1, angle2 in zip(top1.angles, top2.angles):
                if not are_equivalent_connections(angle1, angle2):
                    return False, f"Non equivalent angles {angle1, angle2}"

            for dihedral1, dihedral2 in zip(top1.dihedrals, top2.dihedrals):
                if not are_equivalent_connections(dihedral1, dihedral2):
                    return (
                        False,
                        f"Non equivalent dihedrals, {dihedral1, dihedral2}",
                    )

            for improper1, improper2 in zip(top1.impropers, top2.impropers):
                if not are_equivalent_connections(improper1, improper2):
                    return (
                        False,
                        f"Non equivalent impropers, {improper1, improper2}",
                    )

            for pp_type1, pp_type2 in zip(
                top1.pairpotential_types, top2.pairpotential_types
            ):
                if pp_type1 != pp_type2:
                    return False, f"Pair-PotentialTypes mismatch"

            return True, f"{top1} and {top2} are equivalent"

        return test_topology_equivalence

    @pytest.fixture(scope="session")
    def pairpotentialtype_top(self):
        top = Topology()
        atype1 = AtomType(name="a1")
        atype1.expression = "sigma + epsilon*r"
        atype2 = AtomType(name="a2")
        atype2.expression = "sigma * epsilon * r"
        atom1 = Atom(name="a", atom_type=atype1)
        atom2 = Atom(name="b", atom_type=atype2)
        top.add_site(atom1)
        top.add_site(atom2)
        top.update_topology()

        pptype12 = PairPotentialType(
            name="pp12",
            expression="r + 1",
            independent_variables="r",
            parameters={},
            member_types=tuple(["a1", "a2"]),
        )

        top.add_pairpotentialtype(pptype12)
        return top

    @pytest.fixture(scope="session")
    def residue_top(self):
        top = Topology()
        for i in range(1, 26):
            atom = Atom(
                name=f"atom_{i + 1}",
                residue_number=i % 5,
                residue_name="MY_RES_EVEN" if i % 2 == 0 else f"MY_RES_ODD",
            )
            top.add_site(atom, update_types=False)
        top.update_topology()

        return top
