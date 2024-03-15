import random
from collections import Counter
from operator import attrgetter, itemgetter

import foyer
import mbuild as mb
import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

from gmso.core.views import PotentialFilters
from gmso.external.convert_parmed import from_parmed, to_parmed
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn, has_parmed, import_

pfilter = PotentialFilters.UNIQUE_SORTED_NAMES

if has_parmed:
    pmd = import_("parmed")


@pytest.mark.skipif(not has_parmed, reason="ParmEd is not installed")
class TestConvertParmEd(BaseTest):
    def test_from_parmed_basic(self, angles):
        struc = pmd.load_file(get_fn("ethane.mol2"), structure=True)
        top = from_parmed(struc, refer_type=False)
        for site in top.sites:
            assert site.atom_type is None
        for connection in top.connections:
            assert connection.connection_type is None
        assert top.n_sites == 8
        assert top.n_bonds == 7

        assert top.box is not None
        lengths = u.nm * [0.714, 0.7938, 0.6646]
        assert_allclose_units(top.box.lengths, lengths, rtol=1e-5, atol=1e-8)
        assert_allclose_units(top.box.angles, angles, rtol=1e-5, atol=1e-8)

    def test_from_parmed_parametrized_structure(self, angles):
        struc = pmd.load_file(get_fn("ethane.top"), xyz=get_fn("ethane.gro"))
        top = from_parmed(struc)
        assert top.n_sites == 8
        assert top.n_bonds == 7
        assert top.n_angles == 12
        assert top.n_dihedrals == 9
        assert top.n_connections == 28

        for site in top.sites:
            assert site.atom_type is not None
            assert site.charge is not None

        for connection in top.connections:
            assert connection.connection_type is not None

        assert top.box is not None
        lengths = u.nm * [0.714, 0.7938, 0.6646]
        assert_allclose_units(top.box.lengths, lengths, rtol=1e-5, atol=1e-8)
        assert_allclose_units(top.box.angles, angles, rtol=1e-5, atol=1e-8)

    def test_to_parmed_simple(self):
        struc = pmd.load_file(get_fn("ethane.top"), xyz=get_fn("ethane.gro"))
        struc.title = "Ethane"
        top = from_parmed(struc)
        assert top.name == "Ethane"
        struc_from_top = to_parmed(top, refer_type=False)
        assert struc_from_top.title == "Ethane"

        assert len(struc.atoms) == len(struc_from_top.atoms)
        assert len(struc.bonds) == len(struc_from_top.bonds)
        assert len(struc.angles) == len(struc_from_top.angles)
        assert len(struc.dihedrals) == len(struc_from_top.dihedrals)
        assert len(struc.rb_torsions) == len(struc_from_top.rb_torsions)

    def test_to_parmed_full(self):
        struc = pmd.load_file(get_fn("ethane.top"), xyz=get_fn("ethane.gro"))
        top = from_parmed(struc)
        struc_from_top = to_parmed(top)

        assert struc.bond_types == struc_from_top.bond_types
        assert struc.angle_types == struc_from_top.angle_types
        assert struc.dihedral_types == struc_from_top.dihedral_types
        assert struc.rb_torsion_types == struc_from_top.rb_torsion_types

        # Detail comparisions
        for i in range(len(struc.atoms)):
            assert struc_from_top.atoms[i].name == struc.atoms[i].name
            assert struc_from_top.atoms[i].atom_type == struc.atoms[i].atom_type

        for i in range(len(struc.bonds)):
            assert (
                struc_from_top.bonds[i].atom1.name == struc.bonds[i].atom1.name
                or struc_from_top.bonds[i].atom2.name
                == struc.bonds[i].atom1.name
            )
            assert (
                struc_from_top.bonds[i].atom2.name == struc.bonds[i].atom2.name
                or struc_from_top.bonds[i].atom2.name
                == struc.bonds[i].atom1.name
            )
            assert struc_from_top.bonds[i].type == struc.bonds[i].type

        for i in range(len(struc.angles)):
            assert (
                struc_from_top.angles[i].atom1.name
                == struc.angles[i].atom1.name
            )
            assert (
                struc_from_top.angles[i].atom2.name
                == struc.angles[i].atom2.name
            )
            assert (
                struc_from_top.angles[i].atom3.name
                == struc.angles[i].atom3.name
            )
            assert struc_from_top.angles[i].type == struc.angles[i].type

        for i in range(len(struc.dihedrals)):
            assert (
                struc_from_top.dihedrals[i].atom1.name
                == struc.dihedrals[i].atom1.name
            )
            assert (
                struc_from_top.dihedrals[i].atom2.name
                == struc.dihedrals[i].atom2.name
            )
            assert (
                struc_from_top.dihedrals[i].atom3.name
                == struc.dihedrals[i].atom3.name
            )
            assert (
                struc_from_top.dihedrals[i].atom4.name
                == struc.dihedrals[i].atom4.name
            )
            assert struc_from_top.dihedrals[i].type == struc.dihedrals[i].type

        for i in range(len(struc.rb_torsions)):
            assert (
                struc_from_top.rb_torsions[i].atom1.name
                == struc.rb_torsions[i].atom1.name
                or struc_from_top.rb_torsions[i].atom4.name
                == struc.rb_torsions[i].atom1.name
            )
            assert (
                struc_from_top.rb_torsions[i].atom2.name
                == struc.rb_torsions[i].atom2.name
                or struc_from_top.rb_torsions[i].atom3.name
                == struc.rb_torsions[i].atom2.name
            )
            assert (
                struc_from_top.rb_torsions[i].atom3.name
                == struc.rb_torsions[i].atom3.name
                or struc_from_top.rb_torsions[i].atom2.name
                == struc.rb_torsions[i].atom3.name
            )
            assert (
                struc_from_top.rb_torsions[i].atom4.name
                == struc.rb_torsions[i].atom4.name
                or struc_from_top.rb_torsions[i].atom1.name
                == struc.rb_torsions[i].atom4.name
            )
            assert (
                struc_from_top.rb_torsions[i].type == struc.rb_torsions[i].type
            )

    def test_to_parmed_incompatible_expression(self):
        struc = pmd.load_file(get_fn("ethane.top"), xyz=get_fn("ethane.gro"))
        top = from_parmed(struc)

        with pytest.raises(Exception):
            top.atom_types[0] = "sigma + epsilon"
            struc_from_top = to_parmed(top)

        with pytest.raises(Exception):
            top.bond_types[0] = "k * r_eq"
            struc_from_top = to_parmed(top)

        with pytest.raises(Exception):
            top.angle_types[0] = "k - theta_eq"
            struc_from_top = to_parmed(top)

        with pytest.raises(Exception):
            top.dihedral_types[0] = "c0 - c1 + c2 - c3 + c4 - c5"
            struc_from_top = to_parmed(top)

    def test_to_parmed_loop(
        self, parmed_methylnitroaniline, parmed_chloroethanol, parmed_ethane
    ):
        for struc in [parmed_methylnitroaniline, parmed_chloroethanol]:
            top_from_struc = from_parmed(struc)

            struc_from_top = to_parmed(top_from_struc)

            assert set(struc.bond_types) == set(struc_from_top.bond_types)
            assert set(struc.angle_types) == set(struc_from_top.angle_types)
            assert set(struc.dihedral_types) == set(
                struc_from_top.dihedral_types
            )
            assert set(struc.rb_torsion_types) == set(
                struc_from_top.rb_torsion_types
            )

            # Detail comparisions
            for i in range(len(struc.atoms)):
                assert struc_from_top.atoms[i].name == struc.atoms[i].name
                assert (
                    struc_from_top.atoms[i].atom_type
                    == struc.atoms[i].atom_type
                )

            for i in range(len(struc.bonds)):
                assert (
                    struc_from_top.bonds[i].atom1.name
                    == struc.bonds[i].atom1.name
                )
                assert (
                    struc_from_top.bonds[i].atom2.name
                    == struc.bonds[i].atom2.name
                )
                assert struc_from_top.bonds[i].type == struc.bonds[i].type

            for i in range(len(struc.angles)):
                assert (
                    struc_from_top.angles[i].atom1.name
                    == struc.angles[i].atom1.name
                )
                assert (
                    struc_from_top.angles[i].atom2.name
                    == struc.angles[i].atom2.name
                )
                assert (
                    struc_from_top.angles[i].atom3.name
                    == struc.angles[i].atom3.name
                )
                assert struc_from_top.angles[i].type == struc.angles[i].type

            for i in range(len(struc.dihedrals)):
                assert (
                    struc_from_top.dihedrals[i].atom1.name
                    == struc.dihedrals[i].atom1.name
                )
                assert (
                    struc_from_top.dihedrals[i].atom2.name
                    == struc.dihedrals[i].atom2.name
                )
                assert (
                    struc_from_top.dihedrals[i].atom3.name
                    == struc.dihedrals[i].atom3.name
                )
                assert (
                    struc_from_top.dihedrals[i].atom4.name
                    == struc.dihedrals[i].atom4.name
                )
                assert (
                    struc_from_top.dihedrals[i].type == struc.dihedrals[i].type
                )

            for i in range(len(struc.rb_torsions)):
                assert (
                    struc_from_top.rb_torsions[i].atom1.name
                    == struc.rb_torsions[i].atom1.name
                    or struc_from_top.rb_torsions[i].atom4.name
                    == struc.rb_torsions[i].atom1.name
                )
                assert (
                    struc_from_top.rb_torsions[i].atom2.name
                    == struc.rb_torsions[i].atom2.name
                    or struc_from_top.rb_torsions[i].atom3.name
                    == struc.rb_torsions[i].atom2.name
                )
                assert (
                    struc_from_top.rb_torsions[i].atom3.name
                    == struc.rb_torsions[i].atom3.name
                    or struc_from_top.rb_torsions[i].atom2.name
                    == struc.rb_torsions[i].atom3.name
                )
                assert (
                    struc_from_top.rb_torsions[i].atom4.name
                    == struc.rb_torsions[i].atom4.name
                    or struc_from_top.rb_torsions[i].atom1.name
                    == struc.rb_torsions[i].atom4.name
                )
                assert (
                    struc_from_top.rb_torsions[i].type
                    == struc.rb_torsions[i].type
                )

    def test_residues_info(self, parmed_hexane_box):
        struc = parmed_hexane_box

        top_from_struc = from_parmed(struc)
        assert len(
            top_from_struc.unique_site_labels("residue", name_only=False)
        ) == len(struc.residues)

        for site in top_from_struc.sites:
            assert site.residue[0] == "HEX"
            assert site.residue[1] in list(range(0, 6))

        struc_from_top = to_parmed(top_from_struc)
        assert len(struc_from_top.residues) == len(struc.residues)

        for residue_og, residue_cp in zip(
            struc.residues, struc_from_top.residues
        ):
            assert residue_og.name == residue_cp.name
            assert residue_og.number == residue_cp.number
            assert len(residue_og.atoms) == len(residue_cp.atoms)

    def test_default_residue_info(selfself, parmed_hexane_box):
        struc = parmed_hexane_box
        top_from_struc = from_parmed(struc)
        assert len(
            top_from_struc.unique_site_labels("residue", name_only=False)
        ) == len(struc.residues)

        for site in top_from_struc.sites:
            site.residue = None

        struc_from_top = to_parmed(top_from_struc)
        assert len(struc_from_top.residues) == 1
        assert struc_from_top.residues[0].name == "RES"
        assert len(struc_from_top.atoms) == len(struc.atoms)

    def test_box_info(self, parmed_hexane_box):
        struc = parmed_hexane_box

        top_from_struc = from_parmed(struc)
        assert_allclose_units(
            top_from_struc.box.lengths.to("nm").value,
            [6.0, 6.0, 6.0],
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            top_from_struc.box.angles.to("degree").value,
            [90.0, 90.0, 90.0],
            rtol=1e-5,
            atol=1e-8,
        )

        struc_from_top = to_parmed(top_from_struc)
        assert_allclose_units(
            struc_from_top.box, [60, 60, 60, 90, 90, 90], rtol=1e-5, atol=1e-8
        )

    def test_from_parmed_member_types(self):
        struc = pmd.load_file(get_fn("ethane.top"), xyz=get_fn("ethane.gro"))
        top = from_parmed(struc)
        for potential_types in [
            getattr(top, attr)
            for attr in ["bond_types", "angle_types", "dihedral_types"]
        ]:
            for potential in potential_types:
                assert potential.member_types

    def test_parmed_element(self):
        struc = pmd.load_file(get_fn("ethane.top"), xyz=get_fn("ethane.gro"))
        top = from_parmed(struc)
        for gmso_atom, pmd_atom in zip(top.sites, struc.atoms):
            assert gmso_atom.element.atomic_number == pmd_atom.element

    def test_parmed_element_non_atomistic(self, pentane_ua_parmed):
        top = from_parmed(pentane_ua_parmed, refer_type=False)
        for gmso_atom, pmd_atom in zip(top.sites, pentane_ua_parmed.atoms):
            assert gmso_atom.element is None
            assert pmd_atom.element == 0

    def test_from_parmed_impropers(self):
        mol = "NN-dimethylformamide"
        pmd_structure = pmd.load_file(
            get_fn("{}.top".format(mol)),
            xyz=get_fn("{}.gro".format(mol)),
            parametrize=False,
        )
        assert all(dihedral.improper for dihedral in pmd_structure.dihedrals)
        assert len(pmd_structure.rb_torsions) == 16

        gmso_top = from_parmed(pmd_structure, refer_type=False)
        assert len(gmso_top.impropers) == 2
        for gmso_improper, pmd_improper in zip(
            gmso_top.impropers, pmd_structure.dihedrals
        ):
            pmd_member_names = list(
                atom.name
                for atom in [
                    getattr(pmd_improper, f"atom{j+1}") for j in range(4)
                ]
            )
            gmso_member_names = list(
                map(lambda a: a.name, gmso_improper.connection_members)
            )
            assert pmd_member_names[0] == gmso_member_names[0] and set(
                pmd_member_names[1:]
            ) == set(gmso_member_names[1:])

        assert all(dihedral.improper for dihedral in pmd_structure.dihedrals)
        assert len(pmd_structure.rb_torsions) == 16

    def test_simple_pmd_dihedrals_no_types(self):
        struct = pmd.Structure()
        all_atoms = []
        for j in range(25):
            atom = pmd.Atom(
                atomic_number=j + 1,
                type=f"atom_type_{j + 1}",
                charge=random.randint(1, 10),
                mass=1.0,
            )
            atom.xx, atom.xy, atom.xz = (
                random.random(),
                random.random(),
                random.random(),
            )
            all_atoms.append(atom)
            struct.add_atom(atom, "RES", 1)

        for j in range(10):
            dih = pmd.Dihedral(
                *random.sample(struct.atoms, 4),
                improper=True if j % 2 == 0 else False,
            )
            struct.dihedrals.append(dih)
        gmso_top = from_parmed(struct, refer_type=False)
        assert len(gmso_top.impropers) == 5
        assert len(gmso_top.dihedrals) == 5
        assert len(gmso_top.improper_types) == 0
        assert len(gmso_top.dihedral_types) == 0

    def test_simple_pmd_dihedrals_impropers(self):
        struct = pmd.Structure()
        all_atoms = []
        for j in range(25):
            atype = pmd.AtomType(f"atom_type_{j + 1}", j, 1, atomic_number=12)
            atype.epsilon = atype.rmin = 1
            atom = pmd.Atom(
                atomic_number=12,
                charge=random.randint(1, 10),
                mass=1.0,
            )
            atom.atom_type = atype
            atom.xx, atom.xy, atom.xz = (
                random.random(),
                random.random(),
                random.random(),
            )
            all_atoms.append(atom)
            struct.add_atom(atom, "RES", 1)

        for j in range(10):
            dih = pmd.Dihedral(
                *random.sample(struct.atoms, 4),
                improper=True if j % 2 == 0 else False,
            )
            struct.dihedrals.append(dih)
            dtype = pmd.DihedralType(
                random.random(), random.random(), random.random()
            )
            dih.type = dtype
            struct.dihedral_types.append(dtype)

        gmso_top = from_parmed(struct)
        assert len(gmso_top.impropers) == 5
        assert len(gmso_top.dihedrals) == 5
        assert len(gmso_top.improper_types) == 5
        assert len(gmso_top.dihedral_types) == 5

    def test_pmd_improper_types(self):
        struct = pmd.Structure()
        all_atoms = []
        for j in range(25):
            atype = pmd.AtomType(f"atom_type_{j + 1}", j, 1, atomic_number=12)
            atype.epsilon = atype.rmin = 1
            atom = pmd.Atom(
                atomic_number=j + 1,
                type=f"atom_type_{j + 1}",
                charge=random.randint(1, 10),
                mass=1.0,
            )
            atom.atom_type = atype
            atom.xx, atom.xy, atom.xz = (
                random.random(),
                random.random(),
                random.random(),
            )
            all_atoms.append(atom)
            struct.add_atom(atom, "RES", 1)

        for j in range(10):
            struct.impropers.append(
                pmd.Improper(*random.sample(struct.atoms, 4))
            )
        for improp in struct.impropers:
            improp.type = pmd.ImproperType(random.random(), random.random())
            struct.improper_types.append(improp.type)

        gmso_top = from_parmed(struct)
        assert len(gmso_top.impropers) == len(struct.impropers)
        assert len(gmso_top.improper_types) == len(struct.improper_types)

    def test_pmd_improper_no_types(self):
        struct = pmd.Structure()
        all_atoms = []
        for j in range(25):
            atype = pmd.AtomType(f"atom_type_{j + 1}", j, 1, atomic_number=12)
            atype.epsilon = atype.rmin = 1
            atom = pmd.Atom(
                atomic_number=j + 1,
                type=f"atom_type_{j + 1}",
                charge=random.randint(1, 10),
                mass=1.0,
            )
            atom.atom_atype = atype
            atom.xx, atom.xy, atom.xz = (
                random.random(),
                random.random(),
                random.random(),
            )
            all_atoms.append(atom)
            struct.add_atom(atom, "RES", 1)

        for j in range(10):
            struct.impropers.append(
                pmd.Improper(*random.sample(struct.atoms, 4))
            )

        gmso_top = from_parmed(struct)
        assert len(gmso_top.impropers) == 10
        assert len(gmso_top.improper_types) == 0

    def test_pmd_complex_typed(self, parmed_methylnitroaniline):
        struc = parmed_methylnitroaniline
        top = from_parmed(struc)
        # check connections
        assert top.n_sites == len(struc.atoms)
        assert top.n_bonds == len(struc.bonds)
        assert top.n_angles == len(struc.angles)
        assert top.n_dihedrals == len(struc.rb_torsions)

        # check typing
        assert len(top.atom_types(filter_by=pfilter)) == len(
            Counter(map(attrgetter("atom_type.name"), struc.atoms))
        )
        bonds_list = list(
            map(attrgetter("atom1.type", "atom2.type"), struc.bonds)
        )
        assert len(top.bond_types(filter_by=pfilter)) == len(
            Counter(tuple(sorted(t)) for t in bonds_list)
        )

        angles_list = list(
            map(
                attrgetter("atom1.type", "atom2.type", "atom3.type"),
                struc.angles,
            )
        )
        assert len(top.angle_types(filter_by=pfilter)) == len(
            Counter(
                (t[1], tuple(sorted(itemgetter(*[0, 2])(t))))
                for t in angles_list
            )
        )
        dihedrals_list = list(
            map(
                attrgetter(
                    "atom1.type", "atom2.type", "atom3.type", "atom4.type"
                ),
                struc.rb_torsions,
            )
        )
        # return true if reversal is necessary, false if keep the order
        # order should be from smallest to largest id
        # reverse dihedral order if 1 > 2, or 1=2 and 0>4
        rev_order = lambda x: x[1] > x[2] or (x[1] == x[2] and x[0] > x[3])
        assert len(top.dihedral_types(filter_by=pfilter)) == len(
            Counter(
                tuple(reversed(t)) if rev_order(t) else t
                for t in dihedrals_list
            )
        )

    def test_pmd_complex_impropers_dihedrals(self, parmed_benzene):
        struc = parmed_benzene
        top = from_parmed(struc)
        # check connections
        assert top.n_sites == len(struc.atoms)
        assert top.n_bonds == len(struc.bonds)
        assert top.n_angles == len(struc.angles)
        assert top.n_dihedrals == len(
            [dihedral for dihedral in struc.dihedrals if not dihedral.improper]
        )
        assert top.n_impropers == len(
            [dihedral for dihedral in struc.dihedrals if dihedral.improper]
        )
        # check typing
        assert len(top.atom_types(filter_by=pfilter)) == len(
            Counter(map(attrgetter("atom_type.name"), struc.atoms))
        )
        bonds_list = list(
            map(attrgetter("atom1.type", "atom2.type"), struc.bonds)
        )
        assert len(top.bond_types(filter_by=pfilter)) == len(
            Counter(tuple(sorted(t)) for t in bonds_list)
        )
        angles_list = list(
            map(
                attrgetter("atom1.type", "atom2.type", "atom3.type"),
                struc.angles,
            )
        )
        assert len(top.angle_types(filter_by=pfilter)) == len(
            Counter(
                (t[1], tuple(sorted(itemgetter(*[0, 2])(t))))
                for t in angles_list
            )
        )
        dihedrals = [
            dihedral for dihedral in struc.dihedrals if not dihedral.improper
        ]
        dihedrals_list = list(
            map(
                attrgetter(
                    "atom1.type", "atom2.type", "atom3.type", "atom4.type"
                ),
                dihedrals,
            )
        )
        # return true if reversal is necessary, false if keep the order
        # order should be from smallest to largest id
        # reverse dihedral order if 1 > 2, or 1=2 and 0>4
        rev_order = lambda x: x[1] > x[2] or (x[1] == x[2] and x[0] > x[3])
        assert len(top.dihedral_types(filter_by=pfilter)) == len(
            Counter(
                tuple(reversed(t)) if rev_order(t) else t
                for t in dihedrals_list
            )
        )
        dihedrals = [
            dihedral for dihedral in struc.dihedrals if dihedral.improper
        ]
        impropers_list = list(
            map(
                attrgetter(
                    "atom1.type", "atom2.type", "atom3.type", "atom4.type"
                ),
                dihedrals,
            )
        )
        assert len(top.improper_types(filter_by=pfilter)) == len(
            Counter(t for t in impropers_list)
        )

    def test_pmd_complex_ureybradleys(self, parmed_methylnitroaniline):
        system = mb.Compound()
        first = mb.Particle(name="_CTL2", pos=[-1, 0, 0])
        second = mb.Particle(name="_CL", pos=[0, 0, 0])
        third = mb.Particle(name="_OBL", pos=[1, 0, 0])
        fourth = mb.Particle(name="_OHL", pos=[0, 1, 0])

        system.add([first, second, third, fourth])

        system.add_bond((first, second))
        system.add_bond((second, third))
        system.add_bond((second, fourth))

        ff = foyer.Forcefield(forcefield_files=[get_fn("charmm36_cooh.xml")])
        struc = ff.apply(
            system,
            assert_angle_params=False,
            assert_dihedral_params=False,
            assert_improper_params=False,
        )
        assert len(struc.angles) == 3
        assert len(struc.urey_bradleys) == 2
        top = from_parmed(struc)

        assert top.n_sites == len(struc.atoms)
        assert top.n_bonds == len(struc.bonds)
        assert top.n_angles == len(struc.angles)
        assert top.n_dihedrals == len(struc.rb_torsions)
        # check typing
        assert len(top.atom_types) == len(
            Counter(map(attrgetter("atom_type.name"), struc.atoms))
        )
        bonds_list = list(
            map(attrgetter("atom1.type", "atom2.type"), struc.bonds)
        )
        assert len(top.bond_types(filter_by=pfilter)) == len(
            Counter(tuple(sorted(t)) for t in bonds_list)
        )

        angles_list = list(
            map(
                attrgetter("atom1.type", "atom2.type", "atom3.type"),
                struc.angles,
            )
        )
        assert len(top.angle_types(filter_by=pfilter)) == len(
            Counter(
                (t[1], tuple(sorted(itemgetter(*[0, 2])(t))))
                for t in angles_list
            )
        )
        dihedrals_list = list(
            map(
                attrgetter(
                    "atom1.type", "atom2.type", "atom3.type", "atom4.type"
                ),
                struc.rb_torsions,
            )
        )
        rev_order = lambda x: x[1] > x[2] or (x[1] == x[2] and x[0] > x[3])
        assert len(top.dihedral_types(filter_by=pfilter)) == len(
            Counter(
                tuple(reversed(t)) if rev_order(t) else t
                for t in dihedrals_list
            )
        )

    def test_to_and_from_parmed_with_topology(self, typed_ethane):
        top = typed_ethane
        struc = to_parmed(top)
        top_from_struc = from_parmed(struc)
        assert top.n_sites == top_from_struc.n_sites
        assert top.n_bonds == top_from_struc.n_bonds
        assert top.n_angles == top_from_struc.n_angles
        assert top.n_dihedrals == top_from_struc.n_dihedrals
        assert len(top.atom_types) == len(top_from_struc.atom_types)

    def test_to_and_from_parmed_with_impropers_dihedrals(
        self, typed_methylnitroaniline
    ):
        top = typed_methylnitroaniline
        struc = to_parmed(top)
        top_from_struc = from_parmed(struc)
        assert top.n_sites == top_from_struc.n_sites
        assert top.n_bonds == top_from_struc.n_bonds
        assert top.n_angles == top_from_struc.n_angles
        assert top.n_dihedrals == top_from_struc.n_dihedrals
        assert top.n_impropers == top_from_struc.n_impropers
        assert len(top.atom_types) == len(top_from_struc.atom_types)
        assert len(top.dihedral_types(filter_by=pfilter)) == len(
            top_from_struc.dihedral_types(filter_by=pfilter)
        )
