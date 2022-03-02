import foyer
import mbuild as mb
import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

from gmso.external.convert_parmed import from_parmed, to_parmed
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn, has_parmed, import_

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
            )
            assert (
                struc_from_top.bonds[i].atom2.name == struc.bonds[i].atom2.name
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
            )
            assert (
                struc_from_top.rb_torsions[i].atom2.name
                == struc.rb_torsions[i].atom2.name
            )
            assert (
                struc_from_top.rb_torsions[i].atom3.name
                == struc.rb_torsions[i].atom3.name
            )
            assert (
                struc_from_top.rb_torsions[i].atom4.name
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
                )
                assert (
                    struc_from_top.rb_torsions[i].atom2.name
                    == struc.rb_torsions[i].atom2.name
                )
                assert (
                    struc_from_top.rb_torsions[i].atom3.name
                    == struc.rb_torsions[i].atom3.name
                )
                assert (
                    struc_from_top.rb_torsions[i].atom4.name
                    == struc.rb_torsions[i].atom4.name
                )
                assert (
                    struc_from_top.rb_torsions[i].type
                    == struc.rb_torsions[i].type
                )

    def test_residues_info(self, parmed_hexane_box):
        struc = parmed_hexane_box

        top_from_struc = from_parmed(struc)
        assert len(top_from_struc.subtops) == len(struc.residues)

        for site in top_from_struc.sites:
            assert site.residue_name == "HEX"
            assert site.residue_number in list(range(6))

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
        assert len(top_from_struc.subtops) == len(struc.residues)

        for site in top_from_struc.sites:
            site.residue_name = None
            site.residue_number = None

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
