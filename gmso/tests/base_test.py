import foyer
import mbuild as mb
import numpy as np
import pytest
import unyt as u
from foyer.tests.utils import get_fn

from gmso.core.angle import Angle
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.box import Box
from gmso.core.dihedral import Dihedral
from gmso.core.forcefield import ForceField
from gmso.core.improper import Improper
from gmso.core.pairpotential_type import PairPotentialType
from gmso.core.topology import Topology
from gmso.external import from_mbuild, from_parmed
from gmso.parameterization import apply
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
        return 1 * u.gram / u.mol

    @pytest.fixture
    def box(self):
        return Box(lengths=u.nm * np.ones(3))

    @pytest.fixture
    def top(self):
        return Topology(name="mytop")

    @pytest.fixture
    def benzene_ua(self):
        compound = mb.load(get_fn("benzene_ua.mol2"))
        compound.children[0].name = "BenzeneUA"
        top = from_mbuild(compound)
        top.identify_connections()
        return top

    @pytest.fixture
    def benzene_ua_box(self):
        compound = mb.load(get_fn("benzene_ua.mol2"))
        compound.children[0].name = "BenzeneUA"
        compound_box = mb.packing.fill_box(
            compound=compound, n_compounds=5, density=1
        )
        top = from_mbuild(compound_box)
        top.identify_connections()
        return top

    @pytest.fixture
    def typed_benzene_ua_system(self, benzene_ua_box):
        top = benzene_ua_box
        trappe_benzene = ForceField(get_path("benzene_trappe-ua.xml"))
        top = apply(top=top, forcefields=trappe_benzene, remove_untyped=True)
        return top

    @pytest.fixture
    def benzene_aa(self):
        compound = mb.load(get_fn("benzene.mol2"))
        compound.children[0].name = "BenzeneAA"
        top = from_mbuild(compound)
        top.identify_connections()
        return top

    @pytest.fixture
    def benzene_aa_box(self):
        compound = mb.load(get_fn("benzene.mol2"))
        compound.children[0].name = "BenzeneAA"
        compound_box = mb.packing.fill_box(
            compound=compound, n_compounds=5, density=1
        )
        top = from_mbuild(compound_box)
        top.identify_connections()
        return top

    @pytest.fixture
    def typed_benzene_aa_system(self, benzene_aa_box):
        top = benzene_aa_box
        oplsaa = ForceField("oplsaa")
        top = apply(top=top, forcefields=oplsaa, remove_untyped=True)
        return top

    @pytest.fixture
    def ar_system(self, n_ar_system):
        return from_mbuild(n_ar_system(), parse_label=True)

    @pytest.fixture
    def n_ar_system(self):
        def _topology(n_sites=100):
            ar = mb.Compound(name="Ar", element="Ar")

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
    def spce_water(self):
        spce_comp = mb.lib.molecules.water.WaterSPC()
        spce_ff = ForceField(get_fn("gmso_xmls/test_ffstyles/spce.xml"))
        spce_top = spce_comp.to_gmso()
        spce_top.identify_connections()

        spce_top = apply(spce_top, spce_ff, remove_untyped=True)

        for site in spce_top.sites:
            site.restraint = {
                "kx": 1000 * u.Unit("kJ/(mol*nm**2)"),
                "ky": 1000 * u.Unit("kJ/(mol*nm**2)"),
                "kz": 1000 * u.Unit("kJ/(mol*nm**2)"),
            }

        return spce_top

    @pytest.fixture
    def water_system(self):
        water = Topology(name="water")
        water = water.load(get_path("tip3p.mol2"))
        return water

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
        pmd_ethane = mb_ethane.to_parmed(infer_residues=True)
        pmd_ethane = oplsaa.apply(pmd_ethane)
        top = from_parmed(pmd_ethane)
        top.name = "ethane"
        return top

    @pytest.fixture
    def typed_ethane_opls(self, typed_ethane):
        for dihedral in typed_ethane.dihedrals:
            dihedral.dihedral_type.name = "RyckaertBellemansTorsionPotential"
        top = typed_ethane.convert_potential_styles(
            {"dihedrals": "FourierTorsionPotential"}
        )
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
    def typed_methaneUA(self):
        compound = mb.Compound(name="_CH4", charge=0.0)
        trappe = foyer.Forcefield(name="trappe-ua")
        pmd_structure = trappe.apply(compound)
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
        top.identify_connections()
        ff = ForceField(get_path("tip3p.xml"))
        top = apply(top, ff)
        return top

    @pytest.fixture
    def typed_tip3p_rigid_system(self, water_system):
        top = water_system
        top.identify_connections()
        ff = ForceField(get_path("tip3p-rigid.xml"))
        top = apply(top, ff)
        return top

    @pytest.fixture
    def foyer_fullerene(self):
        from foyer.tests.utils import get_fn

        return ForceField(get_fn("fullerene.xml"))

    @pytest.fixture
    def foyer_periodic(self):
        from foyer.tests.utils import get_fn

        return ForceField(get_fn("oplsaa-periodic.xml"))

    @pytest.fixture
    def foyer_urey_bradley(self):
        from foyer.tests.utils import get_fn

        return ForceField(get_fn("charmm36_cooh.xml"))

    @pytest.fixture
    def foyer_rb_torsion(self):
        from foyer.tests.utils import get_fn

        return ForceField(get_fn("refs-multi.xml"))

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
            equal = lambda x1, x2: (
                u.allclose_units(x1, x2)
                if isinstance(x1, u.unyt_array) and isinstance(x2, u.unyt_array)
                else x1 == x2
            )
            for prop in atom1.model_dump(by_alias=True):
                if not equal(
                    atom2.model_dump().get(prop), atom1.model_dump().get(prop)
                ):
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
            if not np.allclose(top1.scaling_factors, top2.scaling_factors):
                return False, "Mismatch in scaling factors"
            for k, v in top1.molecule_scaling_factors.items():
                if k not in top2.scaling_factors:
                    return False, "Mismatch in scaling factors"
                elif not np.allclose(v, top2.molecule_scaling_factors[k]):
                    return False, "Mismatch in scaling factors"
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
    def labeled_top(self):
        top = Topology()
        for i in range(1, 26):
            atom = Atom(
                name=f"atom_{i + 1}",
                residue=("MY_RES_EVEN" if i % 2 == 0 else f"MY_RES_ODD", i % 5),
                molecule=(
                    "MY_MOL_EVEN" if i % 2 == 0 else f"MY_RES_ODD",
                    i % 5,
                ),
                group="MY_GROUP",
            )
            top.add_site(atom, update_types=False)
        top.update_topology()

        return top

    @pytest.fixture(scope="session")
    def pentane_ua_mbuild(self):
        class PentaneUA(mb.Compound):
            """Create a united-atom pentane compound."""

            def __init__(self):
                super(PentaneUA, self).__init__()
                # Calculate the angle between the two ports
                angle = np.deg2rad(114)
                x = 0
                y = 0.077
                z = 0
                vec = [
                    x,
                    y * np.cos(angle) - z * np.sin(angle),
                    y * np.sin(angle) + z * np.cos(angle),
                ]
                # Create the end group compound
                ch3 = mb.Compound()
                ch3.add(mb.Particle(name="_CH3"))
                ch3.add(mb.Port(anchor=ch3[0]), "up")
                ch3["up"].translate([x, y, z])
                # Create the internal monomer
                ch2 = mb.Compound()
                ch2.add(mb.Particle(name="_CH2"))
                ch2.add(mb.Port(anchor=ch2[0]), "up")
                ch2["up"].translate([x, y, z])
                ch2.add(mb.Port(anchor=ch2[0], orientation=vec), "down")
                ch2["down"].translate(vec)
                pentane = mb.recipes.Polymer(
                    monomers=[ch2], end_groups=[ch3, mb.clone(ch3)]
                )
                pentane.build(n=3)
                self.add(pentane, label="PNT")

        return PentaneUA()

    @pytest.fixture(scope="session")
    def pentane_ua_parmed(self, pentane_ua_mbuild):
        return mb.conversion.to_parmed(pentane_ua_mbuild)

    @pytest.fixture(scope="session")
    def pentane_ua_gmso(self, pentane_ua_mbuild):
        return from_mbuild(pentane_ua_mbuild)

    @pytest.fixture(scope="session")
    def hierarchical_compound(self):
        # Build Polymer
        monomer = mb.load("CCO", smiles=True)
        monomer.name = "monomer"
        polymer = mb.lib.recipes.Polymer()
        polymer.add_monomer(monomer, indices=(3, 7))
        polymer.build(n=10)
        polymer.name = "polymer"

        # Build Solvent 1
        cyclopentane = mb.load("C1CCCC1", smiles=True)
        cyclopentane.name = "cyclopentane"

        # Build Solvent 2
        water = mb.load("O", smiles=True)
        water.name = "water"

        # Build Partitioned Box
        filled_box1 = mb.packing.solvate(
            solvent=cyclopentane,
            solute=polymer,
            box=mb.Box([5, 5, 5]),
            n_solvent=5,
        )
        filled_box1.name = "sol1"
        filled_box2 = mb.packing.fill_box(
            compound=water,
            box=mb.Box([5, 5, 5]),
            n_compounds=5,
        )
        filled_box2.name = "sol2"
        partitioned_box = mb.Compound()
        partitioned_box.add(filled_box1)
        partitioned_box.add(filled_box2)
        partitioned_box.name = "Topology"
        return partitioned_box

    @pytest.fixture(scope="session")
    def hierarchical_top(self, hierarchical_compound):
        top = from_mbuild(hierarchical_compound)  # Create GMSO topology
        top.identify_connections()
        return top

    @pytest.fixture
    def ethane_gomc(self):
        ethane_gomc = mb.load("CC", smiles=True)
        ethane_gomc.name = "ETH"

        return ethane_gomc

    @pytest.fixture
    def ethanol_gomc(self):
        ethanol_gomc = mb.load("CCO", smiles=True)
        ethanol_gomc.name = "ETO"

        return ethanol_gomc

    @pytest.fixture
    def methane_ua_gomc(self):
        methane_ua_gomc = mb.Compound(name="_CH4")

        return methane_ua_gomc

    @pytest.fixture
    def parmed_benzene(self):
        untyped_benzene = mb.load(get_fn("benzene.mol2"))
        ff_improper = foyer.Forcefield(
            forcefield_files=get_fn("improper_dihedral.xml")
        )
        benzene = ff_improper.apply(
            untyped_benzene, assert_dihedral_params=False
        )
        return benzene

    @pytest.fixture
    def benzeneTopology(self):
        untyped_benzene = mb.load(get_fn("benzene.mol2"))
        top_benzene = untyped_benzene.to_gmso()
        ff_improper = ForceField(get_fn("benzeneaa_improper.xml"))
        return apply(
            top_benzene,
            ff_improper,
            identify_connections=True,
            ignore_params=[],
        )

    # TODO: now
    # add in some fixtures for (connects), amber

    @pytest.fixture
    def harmonic_parmed_types_charmm(self):
        from mbuild.formats.lammpsdata import write_lammpsdata

        system = mb.Compound()
        first = mb.Particle(name="_CTL2", pos=[-1, 0, 0])
        second = mb.Particle(name="_CL", pos=[0, 0, 0])
        third = mb.Particle(name="_OBL", pos=[1, 0, 0])
        fourth = mb.Particle(name="_OHL", pos=[0, 1, 0])

        system.add([first, second, third, fourth])

        system.add_bond((first, second))
        system.add_bond((second, third))
        system.add_bond((second, fourth))

        ff = foyer.Forcefield(forcefield_files=[get_path("charmm36_cooh.xml")])
        struc = ff.apply(
            system,
            assert_angle_params=False,
            assert_dihedral_params=False,
            assert_improper_params=False,
        )
        return struc

    @pytest.fixture
    def gaff_forcefield(self):
        return ForceField(get_fn("gmso_xmls/test_ffstyles/gaff.xml"))

    @pytest.fixture
    def oplsaa_forcefield(self):
        return ForceField("oplsaa")
