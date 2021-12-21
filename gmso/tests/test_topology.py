from copy import deepcopy

import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType
from gmso.core.box import Box
from gmso.core.dihedral import Dihedral
from gmso.core.dihedral_type import DihedralType
from gmso.core.improper import Improper
from gmso.core.improper_type import ImproperType
from gmso.core.pairpotential_type import PairPotentialType
from gmso.core.subtopology import SubTopology
from gmso.core.topology import Topology
from gmso.exceptions import GMSOError
from gmso.external.convert_parmed import from_parmed
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn, has_parmed, import_

if has_parmed:
    pmd = import_("parmed")


class TestTopology(BaseTest):
    def test_new_topology(self):
        top = Topology(name="mytop")
        assert top.name == "mytop"

    def test_empty_name(self):
        top = Topology(name="")
        assert top.name == "Topology"

    def test_change_comb_rule(self):
        top = Topology()
        assert top.combining_rule == "lorentz"
        top.combining_rule = "geometric"
        assert top.combining_rule == "geometric"
        with pytest.raises(GMSOError):
            top.combining_rule = "kong"

    def test_add_site(self):
        top = Topology()
        site = Atom(name="site")

        assert top.n_sites == 0
        top.add_site(site)
        assert top.n_sites == 1

    def test_add_connection(self):
        top = Topology()
        atom1 = Atom(name="atom1")
        atom2 = Atom(name="atom2")
        connect = Bond(connection_members=[atom1, atom2])

        top.add_connection(connect)
        top.add_site(atom1)
        top.add_site(atom2)

        assert len(top.connections) == 1

    def test_add_box(self):
        top = Topology()
        box = Box(2 * u.nm * np.ones(3))

        assert top.box is None
        top.box = box
        assert top.box is not None
        assert_allclose_units(
            top.box.lengths, u.nm * 2 * np.ones(3), rtol=1e-5, atol=1e-8
        )

    def test_positions_dtype(self):
        top = Topology()
        atom1 = Atom(name="atom1", position=[0.0, 0.0, 0.0])
        top.add_site(atom1)

        assert set([type(site.position) for site in top.sites]) == {
            u.unyt_array
        }
        assert set([site.position.units for site in top.sites]) == {u.nm}

        assert top.positions.dtype == float
        assert top.positions.units == u.nm
        assert isinstance(top.positions, u.unyt_array)

    def test_eq_types(self, top, box):
        assert top != box

        diff_name = deepcopy(top)
        diff_name.name = "othertop"
        assert top != diff_name

    def test_eq_sites(self, top, charge):
        ref = deepcopy(top)
        wrong_n_sites = deepcopy(top)
        assert top != wrong_n_sites
        ref.add_site(Atom())
        assert ref != wrong_n_sites

        ref = deepcopy(top)
        wrong_position = deepcopy(top)
        ref.add_site(Atom(position=u.nm * [0, 0, 0]))
        wrong_position.add_site(Atom(position=u.nm * [1, 1, 1]))
        assert top != wrong_position

        ref = deepcopy(top)
        wrong_charge = deepcopy(top)
        ref.add_site(Atom(charge=charge))
        wrong_charge.add_site(Atom(charge=-1 * charge))
        assert ref != wrong_charge

        ref = deepcopy(top)
        wrong_atom_type = deepcopy(top)
        ref.add_site(Atom(atom_type=AtomType(expression="epsilon*sigma*r")))
        wrong_atom_type.add_site(Atom(atom_type=AtomType(expression="sigma*r")))
        assert ref != wrong_atom_type

    @pytest.mark.skipif(not has_parmed, reason="ParmEd is not installed")
    def test_eq_bonds(self):
        ref = pmd.load_file(get_fn("ethane.top"), xyz=get_fn("ethane.gro"))

        missing_bond = deepcopy(ref)
        missing_bond.bonds[0].delete()

        assert ref != missing_bond

        bad_bond_type = deepcopy(ref)
        bad_bond_type.bond_types[0].k = 22

        assert ref != bad_bond_type

    @pytest.mark.skipif(not has_parmed, reason="ParmEd is not installed")
    def test_eq_angles(self):
        ref = pmd.load_file(get_fn("ethane.top"), xyz=get_fn("ethane.gro"))

        missing_angle = deepcopy(ref)
        missing_angle.angles[0].delete()

        assert ref != missing_angle

        bad_angle_type = deepcopy(ref)
        bad_angle_type.angle_types[0].k = 22

        assert ref != bad_angle_type

    @pytest.mark.skipif(not has_parmed, reason="ParmEd is not installed")
    def test_eq_dihedrals(self):
        ref = pmd.load_file(get_fn("ethane.top"), xyz=get_fn("ethane.gro"))

        missing_dihedral = deepcopy(ref)
        missing_dihedral.rb_torsions[0].delete()

        assert ref != missing_dihedral

        bad_dihedral_type = deepcopy(ref)
        bad_dihedral_type.rb_torsion_types[0].k = 22

        assert ref != bad_dihedral_type

    @pytest.mark.skipif(not has_parmed, reason="ParmEd is not installed")
    def test_eq_overall(self):
        ref = pmd.load_file(get_fn("ethane.top"), xyz=get_fn("ethane.gro"))

        top1 = from_parmed(ref)
        top2 = from_parmed(ref)

        assert top1 != top2

    def test_add_untyped_site_update(self):
        untyped_site = Atom(atom_type=None)

        top = Topology()
        assert len(top.atom_types) == 0
        top.add_site(untyped_site, update_types=False)
        assert len(top.atom_types) == 0

        top = Topology()
        assert len(top.atom_types) == 0
        top.add_site(untyped_site, update_types=True)
        assert len(top.atom_types) == 0

    def test_add_typed_site_update(self):
        typed_site = Atom(atom_type=AtomType())

        top = Topology()
        assert len(top.atom_types) == 0
        top.add_site(typed_site, update_types=False)
        assert len(top.atom_types) == 0

        top = Topology()
        assert len(top.atom_types) == 0
        top.add_site(typed_site, update_types=True)
        assert len(top.atom_types) == 1

    def test_add_untyped_bond_update(self):
        atom1 = Atom(atom_type=None)
        atom2 = Atom(atom_type=None)
        bond = Bond(connection_members=[atom1, atom2], bond_type=None)

        top = Topology()
        assert len(top.bond_types) == 0
        top.add_connection(bond, update_types=False)
        assert len(top.bond_types) == 0

        top = Topology()
        assert len(top.bond_types) == 0
        top.add_connection(bond, update_types=True)
        assert len(top.bond_types) == 0

    def test_add_typed_bond_update(self):
        atom1 = Atom(atom_type=None)
        atom2 = Atom(atom_type=None)
        bond = Bond(connection_members=[atom1, atom2], bond_type=BondType())

        top = Topology()
        top.add_site(atom1)
        top.add_site(atom2)
        top.add_connection(bond, update_types=False)
        assert len(top.connection_types) == 0

        top = Topology()
        top.add_connection(bond, update_types=True)
        assert len(top.bond_types) == 1

    def test_top_update(self):
        top = Topology()
        top.update_topology()
        assert top.n_sites == 0
        assert len(top.atom_types) == 0
        assert len(top.atom_type_expressions) == 0
        assert top.n_connections == 0
        assert len(top.connection_types) == 0
        assert len(top.connection_type_expressions) == 0

        atomtype = AtomType()
        atom1 = Atom(name="atom1", atom_type=atomtype)
        top.add_site(atom1)
        atom2 = Atom(name="atom2", atom_type=atomtype)
        top.add_site(atom2)

        assert top.n_sites == 2
        assert len(top.atom_types) == 1
        assert len(top.atom_type_expressions) == 1
        assert top.n_connections == 0
        assert len(top.connection_types) == 0
        assert len(top.connection_type_expressions) == 0

        ctype = BondType()
        connection_12 = Bond(connection_members=[atom1, atom2], bond_type=ctype)
        top.add_connection(connection_12)

        assert top.n_sites == 2
        assert len(top.atom_types) == 1
        assert len(top.atom_type_expressions) == 1
        assert top.n_connections == 1
        assert len(top.connection_types) == 1
        assert len(top.connection_type_expressions) == 1

        atom1.atom_type = AtomType(expression="sigma*epsilon*r")
        assert top.n_sites == 2
        assert len(top.atom_types) == 1
        assert len(top.atom_type_expressions) == 1
        assert top.n_connections == 1
        assert len(top.connection_types) == 1
        assert len(top.connection_type_expressions) == 1
        top.update_atom_types()
        assert top.n_sites == 2
        assert len(top.atom_types) == 2
        assert len(top.atom_type_expressions) == 2
        assert top.n_connections == 1
        assert len(top.connection_types) == 1
        assert len(top.connection_type_expressions) == 1

    def test_atomtype_update(self):
        top = Topology()

        assert top.n_sites == 0
        assert top.n_bonds == 0
        assert top.n_connections == 0

        atype1 = AtomType(expression="sigma + epsilon*r")
        atype2 = AtomType(expression="sigma * epsilon*r")
        atom1 = Atom(name="a", atom_type=atype1)
        atom2 = Atom(name="b", atom_type=atype2)
        top.add_site(atom1)
        top.add_site(atom2)

        assert top.n_sites == 2
        assert len(top.atom_types) == 2
        assert len(top.atom_type_expressions) == 2

    def test_bond_bondtype_update(self):
        top = Topology()

        atype1 = AtomType(expression="sigma + epsilon*r")
        atype2 = AtomType(expression="sigma * epsilon*r")
        atom1 = Atom(name="a", atom_type=atype1)
        atom2 = Atom(name="b", atom_type=atype2)
        btype = BondType()
        bond = Bond(connection_members=[atom1, atom2], bond_type=btype)
        top.add_site(atom1)
        top.add_site(atom2)
        top.add_connection(bond)

        assert top.n_bonds == 1
        assert len(top.bond_types) == 1
        assert len(top.bond_type_expressions) == 1

    def test_angle_angletype_update(self):
        top = Topology()

        atype1 = AtomType(expression="sigma + epsilon*r")
        atype2 = AtomType(expression="sigma * epsilon*r")
        atom1 = Atom(name="a", atom_type=atype1)
        atom2 = Atom(name="b", atom_type=atype2)
        atom3 = Atom(name="c", atom_type=atype2)
        atype = AngleType()
        angle = Angle(
            connection_members=[atom1, atom2, atom3],
            angle_type=atype,
            name="angle_name",
        )
        top.add_site(atom1)
        top.add_site(atom2)
        top.add_site(atom3)
        top.add_connection(angle)

        assert top.n_angles == 1
        assert len(top.angle_types) == 1
        assert len(top.angle_type_expressions) == 1
        assert len(top.atom_type_expressions) == 2

    def test_dihedral_dihedraltype_update(self):
        top = Topology()

        atype1 = AtomType(expression="sigma + epsilon*r")
        atype2 = AtomType(expression="sigma * epsilon*r")
        atom1 = Atom(name="a", atom_type=atype1)
        atom2 = Atom(name="b", atom_type=atype2)
        atom3 = Atom(name="c", atom_type=atype2)
        atom4 = Atom(name="d", atom_type=atype1)
        atype = DihedralType()
        dihedral = Dihedral(
            connection_members=[atom1, atom2, atom3, atom4], dihedral_type=atype
        )
        top.add_site(atom1)
        top.add_site(atom2)
        top.add_site(atom3)
        top.add_site(atom4)
        top.add_connection(dihedral)

        assert top.n_dihedrals == 1
        assert len(top.dihedral_types) == 1
        assert len(top.dihedral_type_expressions) == 1
        assert len(top.atom_type_expressions) == 2

    def test_improper_impropertype_update(self):
        top = Topology()

        atype1 = AtomType(expression="sigma + epsilon*r")
        atype2 = AtomType(expression="sigma * epsilon*r")
        atom1 = Atom(name="a", atom_type=atype1)
        atom2 = Atom(name="b", atom_type=atype2)
        atom3 = Atom(name="c", atom_type=atype2)
        atom4 = Atom(name="d", atom_type=atype1)
        atype = ImproperType()
        improper = Improper(
            connection_members=[atom1, atom2, atom3, atom4], improper_type=atype
        )
        top.add_site(atom1)
        top.add_site(atom2)
        top.add_site(atom3)
        top.add_site(atom4)
        top.add_connection(improper)

        assert top.n_impropers == 1
        assert len(top.improper_types) == 1
        assert len(top.improper_type_expressions) == 1
        assert len(top.atom_type_expressions) == 2

    def test_pairpotential_pairpotentialtype_update(
        self, pairpotentialtype_top
    ):
        assert len(pairpotentialtype_top.pairpotential_types) == 1
        pptype12 = pairpotentialtype_top.pairpotential_types[0]
        assert pairpotentialtype_top._pairpotential_types_idx[pptype12] == 0

        pairpotentialtype_top.remove_pairpotentialtype(["a1", "a2"])
        assert len(pairpotentialtype_top.pairpotential_types) == 0

    def test_add_subtopology(self):
        top = Topology()
        subtop = SubTopology()

        assert top.n_subtops == 0
        top.add_subtopology(subtop)
        assert top.n_subtops == 1

    def test_parametrization(self):
        top = Topology()

        assert top.typed == False
        top.add_site(Atom(atom_type=AtomType()))

        assert top.typed == True
        assert top.is_typed() == True
        assert top.typed == True

    def test_parametrization_setter(self):
        top = Topology()

        assert top.typed == False
        assert top.is_typed() == False
        top.typed = True
        assert top.typed == True
        assert top.is_typed() == False

    def test_topology_atom_type_changes(self):
        top = Topology()
        for i in range(100):
            site = Atom(name="site{}".format(i))
            atom_type = AtomType(name="atom_type{}".format(i % 10))
            site.atom_type = atom_type
            top.add_site(site, update_types=False)
        top.update_topology()
        assert len(top.atom_types) == 10
        top.sites[0].atom_type.name = "atom_type_changed"
        assert id(top.sites[0].atom_type) == id(top.sites[10].atom_type)
        assert top.sites[10].atom_type.name == "atom_type_changed"
        assert top.is_typed()

    def test_add_duplicate_connected_atom(self):
        top = Topology()
        atom1 = Atom(name="AtomA")
        atom2 = Atom(name="AtomB")
        bond = Bond(connection_members=[atom1, atom2])
        bond_eq = Bond(connection_members=[atom1, atom2])

        top.add_connection(bond)
        top.add_connection(bond_eq)
        top.update_topology()
        assert top.n_connections == 1

    def test_topology_get_index(self):
        top = Topology()
        conn_members = [Atom() for _ in range(10)]
        for atom in conn_members:
            top.add_site(atom)

        for i in range(5):
            top.add_connection(
                Bond(connection_members=[conn_members[i], conn_members[i + 1]])
            )
            top.add_connection(
                Angle(
                    connection_members=[
                        conn_members[i],
                        conn_members[i + 1],
                        conn_members[i + 2],
                    ]
                )
            )
            top.add_connection(
                Dihedral(
                    connection_members=[
                        conn_members[i],
                        conn_members[i + 1],
                        conn_members[i + 2],
                        conn_members[i + 3],
                    ]
                )
            )
            top.add_connection(
                Improper(
                    connection_members=[
                        conn_members[i],
                        conn_members[i + 1],
                        conn_members[i + 2],
                        conn_members[i + 3],
                    ]
                )
            )

        a_atom = Atom()
        a_bond = Bond(connection_members=[conn_members[6], conn_members[7]])
        an_angle = Angle(
            connection_members=[
                conn_members[6],
                conn_members[7],
                conn_members[8],
            ]
        )
        a_dihedral = Dihedral(
            connection_members=[
                conn_members[6],
                conn_members[7],
                conn_members[8],
                conn_members[9],
            ]
        )
        an_improper = Improper(
            connection_members=[
                conn_members[6],
                conn_members[7],
                conn_members[8],
                conn_members[9],
            ]
        )

        top.add_site(a_atom)
        top.add_connection(a_bond)
        top.add_connection(an_angle)
        top.add_connection(a_dihedral)
        top.add_connection(an_improper)

        assert top.get_index(a_atom) == 10
        assert top.get_index(a_bond) == 5
        assert top.get_index(an_angle) == 5
        assert top.get_index(a_dihedral) == 5
        assert top.get_index(an_improper) == 5

    def test_topology_get_index_wrong_member_type(self):
        top = Topology()
        with pytest.raises(TypeError):
            top.get_index(object())

    def test_topology_get_index_non_existing_member(self):
        top = Topology()
        site = Atom()
        with pytest.raises(ValueError):
            top.get_index(site)

    def test_topology_get_index_atom_type(self, typed_water_system):
        assert (
            typed_water_system.get_index(typed_water_system.sites[0].atom_type)
            == 0
        )
        assert (
            typed_water_system.get_index(typed_water_system.sites[1].atom_type)
            == 1
        )

    def test_topology_get_index_atom_type_after_change(
        self, typed_water_system
    ):
        typed_water_system.sites[0].atom_type.name = "atom_type_changed_name"
        assert (
            typed_water_system.get_index(typed_water_system.sites[0].atom_type)
            == 1
        )
        assert (
            typed_water_system.get_index(typed_water_system.sites[1].atom_type)
            == 0
        )

    def test_topology_get_index_bond_type(self, typed_methylnitroaniline):
        assert (
            typed_methylnitroaniline.get_index(
                typed_methylnitroaniline.bonds[0].connection_type
            )
            == 0
        )
        assert isinstance(
            typed_methylnitroaniline.get_index(
                typed_methylnitroaniline.bonds[-1].connection_type
            ),
            int,
        )

    def test_topology_get_index_bond_type_after_change(
        self, typed_methylnitroaniline
    ):
        typed_methylnitroaniline.bonds[0].connection_type.name = "changed name"
        assert (
            typed_methylnitroaniline.get_index(
                typed_methylnitroaniline.bonds[0].connection_type
            )
            != 0
        )

    def test_topology_get_index_angle_type(self, typed_chloroethanol):
        assert (
            typed_chloroethanol.get_index(
                typed_chloroethanol.angles[0].connection_type
            )
            == 0
        )
        assert (
            typed_chloroethanol.get_index(
                typed_chloroethanol.angles[5].connection_type
            )
            == 1
        )

    def test_topology_get_index_angle_type_after_change(
        self, typed_methylnitroaniline
    ):
        angle_type_to_test = typed_methylnitroaniline.angles[0].connection_type
        prev_idx = typed_methylnitroaniline.get_index(angle_type_to_test)
        typed_methylnitroaniline.angles[0].connection_type.name = "changed name"
        assert (
            typed_methylnitroaniline.get_index(angle_type_to_test) != prev_idx
        )

    def test_topology_get_index_dihedral_type(self, typed_chloroethanol):
        assert (
            typed_chloroethanol.get_index(
                typed_chloroethanol.dihedrals[0].connection_type
            )
            == 0
        )
        assert (
            typed_chloroethanol.get_index(
                typed_chloroethanol.dihedrals[5].connection_type
            )
            == 3
        )

    def test_topology_get_index_dihedral_type_after_change(
        self, typed_methylnitroaniline
    ):
        dihedral_type_to_test = typed_methylnitroaniline.dihedrals[
            0
        ].connection_type
        prev_idx = typed_methylnitroaniline.get_index(dihedral_type_to_test)
        typed_methylnitroaniline.dihedrals[
            0
        ].connection_type.name = "changed name"
        assert (
            typed_methylnitroaniline.get_index(dihedral_type_to_test)
            != prev_idx
        )

    def test_topology_get_bonds_for(self, typed_methylnitroaniline):
        site = list(typed_methylnitroaniline.sites)[0]
        converted_bonds_list = typed_methylnitroaniline._get_bonds_for(site)
        top_bonds_containing_site = []
        for bond in typed_methylnitroaniline.bonds:
            if site in bond.connection_members:
                assert bond in converted_bonds_list
                top_bonds_containing_site.append(bond)
        assert len(top_bonds_containing_site) == len(converted_bonds_list)

    def test_topology_get_angles_for(self, typed_methylnitroaniline):
        site = list(typed_methylnitroaniline.sites)[0]
        converted_angles_list = typed_methylnitroaniline._get_angles_for(site)
        top_angles_containing_site = []
        for angle in typed_methylnitroaniline.angles:
            if site in angle.connection_members:
                assert angle in converted_angles_list
                top_angles_containing_site.append(angle)
        assert len(top_angles_containing_site) == len(converted_angles_list)

    def test_topology_get_dihedrals_for(self, typed_methylnitroaniline):
        site = list(typed_methylnitroaniline.sites)[0]
        converted_dihedrals_list = typed_methylnitroaniline._get_dihedrals_for(
            site
        )
        top_dihedrals_containing_site = []
        for dihedral in typed_methylnitroaniline.dihedrals:
            if site in dihedral.connection_members:
                assert dihedral in converted_dihedrals_list
                top_dihedrals_containing_site.append(dihedral)
        assert len(top_dihedrals_containing_site) == len(
            converted_dihedrals_list
        )

    def test_topology_scale_factors(self, typed_methylnitroaniline):
        sf = typed_methylnitroaniline.scaling_factors
        assert np.allclose(sf["nonBonded12Scale"], 0.0)
        assert np.allclose(sf["nonBonded13Scale"], 0.0)
        assert np.allclose(sf["nonBonded14Scale"], 0.5)
        assert np.allclose(sf["electrostatics12Scale"], 0.0)
        assert np.allclose(sf["electrostatics13Scale"], 0.0)
        assert np.allclose(sf["electrostatics14Scale"], 0.5)

    def test_topology_change_scale_factors(self, typed_methylnitroaniline):
        typed_methylnitroaniline.scaling_factors = {
            "nonBonded12Scale": 0.5,
            "nonBonded13Scale": 0.5,
            "nonBonded14Scale": 1.0,
            "electrostatics12Scale": 1.0,
            "electrostatics13Scale": 1.0,
            "electrostatics14Scale": 1.0,
        }
        sf = typed_methylnitroaniline.scaling_factors
        assert np.allclose(sf["nonBonded12Scale"], 0.5)
        assert np.allclose(sf["nonBonded13Scale"], 0.5)
        assert np.allclose(sf["nonBonded14Scale"], 1.0)
        assert np.allclose(sf["electrostatics12Scale"], 1.0)
        assert np.allclose(sf["electrostatics13Scale"], 1.0)
        assert np.allclose(sf["electrostatics14Scale"], 1.0)
        typed_methylnitroaniline.scaling_factors["nonBonded12Scale"] = 1.0
        assert np.allclose(sf["nonBonded12Scale"], 1.0)

    def test_topology_invalid_scaling_factors(self, typed_methylnitroaniline):
        with pytest.raises(GMSOError):
            typed_methylnitroaniline.scaling_factors = (0.5, 1.0)
        with pytest.raises(GMSOError):
            typed_methylnitroaniline.scaling_factors = {"lj_12": 0.0}

    def test_is_typed_check(self, typed_chloroethanol):
        groups = [
            "sites",
            "bonds",
            "angles",
            "dihedrals",
            "impropers",
            "topology",
        ]
        for group in groups:
            assert typed_chloroethanol.is_fully_typed(group=group)

        with pytest.raises(ValueError):
            typed_chloroethanol.is_fully_typed(group="foo")

    def test_cget_untyped(self, typed_chloroethanol):
        # Note impropers list is empty, and hence is not tested here
        groups = ["sites", "bonds", "angles", "dihedrals"]
        clone = deepcopy(typed_chloroethanol)
        clone.sites[0].atom_type = None
        clone.bonds[0].bond_type = None
        clone.angles[0].angle_type = None
        clone.dihedrals[0].dihedral_type = None

        full_opt = clone.get_untyped(group="topology")
        for group in groups:
            assert not clone.is_fully_typed(group=group)
            assert len(clone.get_untyped(group=group)[group]) == 1
            assert len(full_opt[group]) == 1

        # Get multiple untyped
        untyped_sites_angles = clone.get_untyped(group=["sites", "angles"])
        assert len(untyped_sites_angles) == 2
        assert len(untyped_sites_angles["sites"]) == 1
        assert len(untyped_sites_angles["angles"]) == 1

        with pytest.raises(ValueError):
            clone.get_untyped(group="foo")

    def test_iter_sites(self, residue_top):
        for site in residue_top.iter_sites("residue_name", "MY_RES_EVEN"):
            assert site.residue_name == "MY_RES_EVEN"

        for site in residue_top.iter_sites("residue_name", "MY_RES_ODD"):
            assert site.residue_name == "MY_RES_ODD"

        sites = list(residue_top.iter_sites("residue_number", 4))
        assert len(sites) == 5

    def test_iter_sites_non_iterable_attribute(self, residue_top):
        with pytest.raises(ValueError):
            for site in residue_top.iter_sites("atom_type", "abc"):
                pass

    def test_iter_sites_none(self, residue_top):
        with pytest.raises(ValueError):
            for site in residue_top.iter_sites("residue_name", None):
                pass

    def test_iter_sites_by_residue_name(self, pairpotentialtype_top):
        assert (
            len(list(pairpotentialtype_top.iter_sites_by_residue_name("AAA")))
            == 0
        )

    def test_iter_sites_by_residue_number(self, residue_top):
        sites = list(residue_top.iter_sites_by_residue_number(4))
        assert len(sites) == 5
