import unyt as u
import pytest

import mbuild as mb
import foyer

from gmso.external.convert_parmed import from_parmed, to_parmed
from gmso.tests.base_test import BaseTest
from gmso.utils.testing import allclose
from gmso.utils.io import get_fn, import_, has_parmed


if has_parmed:
    pmd = import_('parmed')

@pytest.mark.skipif(not has_parmed, reason="ParmEd is not installed")
class TestConvertParmEd(BaseTest):
    def test_from_parmed_basic(self, angles):
        struc = pmd.load_file(get_fn('ethane.mol2'), structure=True)
        top = from_parmed(struc, refer_type=False)
        for site in top.sites:
            assert site.atom_type is None
        for connection in top.connections:
            assert connection.connection_type is None
        assert top.n_sites == 8
        assert top.n_bonds == 7

        assert top.box is not None
        lengths = u.nm * [0.714, 0.7938, 0.6646]
        assert allclose(top.box.lengths, lengths)
        assert allclose(top.box.angles, angles)

    def test_from_parmed_parametrized_structure(self, angles):
        struc = pmd.load_file(get_fn('ethane.top'), xyz=get_fn('ethane.gro'))
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
        assert allclose(top.box.lengths, lengths)
        assert allclose(top.box.angles, angles)

    def test_to_parmed_simple(self):
        struc = pmd.load_file(get_fn('ethane.top'), xyz=get_fn('ethane.gro'))
        top = from_parmed(struc)
        struc_from_top = to_parmed(top, refer_type=False)

        assert len(struc.atoms) == len(struc_from_top.atoms)
        assert len(struc.bonds) == len(struc_from_top.bonds)
        assert len(struc.angles) == len(struc_from_top.angles)
        assert len(struc.dihedrals) == len(struc.dihedrals)

    def test_to_parmed_full(self):
        struc = pmd.load_file(get_fn('ethane.top'), xyz=get_fn('ethane.gro'))
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
            assert struc_from_top.bonds[i].atom1.name == struc.bonds[i].atom1.name
            assert struc_from_top.bonds[i].atom2.name == struc.bonds[i].atom2.name
            assert struc_from_top.bonds[i].type == struc.bonds[i].type

        for i in range(len(struc.angles)):
            assert struc_from_top.angles[i].atom1.name == struc.angles[i].atom1.name
            assert struc_from_top.angles[i].atom2.name == struc.angles[i].atom2.name
            assert struc_from_top.angles[i].atom3.name == struc.angles[i].atom3.name
            assert struc_from_top.angles[i].type == struc.angles[i].type

        for i in range(len(struc.dihedrals)):
            assert struc_from_top.dihedrals[i].atom1.name == struc.dihedrals[i].atom1.name
            assert struc_from_top.dihedrals[i].atom2.name == struc.dihedrals[i].atom2.name
            assert struc_from_top.dihedrals[i].atom3.name == struc.dihedrals[i].atom3.name
            assert struc_from_top.dihedrals[i].atom4.name == struc.dihedrals[i].atom4.name
            assert struc_from_top.dihedrals[i].type == struc.dihedrals[i].type

    def test_to_parmed_incompatible_expression(self):
        struc = pmd.load_file(get_fn('ethane.top'), xyz=get_fn('ethane.gro'))
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

    def test_to_parmed_loop(self, parmed_methylnitroaniline,
                                  parmed_chloroethanol):
        for struc in [parmed_methylnitroaniline, parmed_chloroethanol]:
            top_from_struc = from_parmed(struc)

            struc_from_top = to_parmed(top_from_struc)

            assert set(struc.bond_types) == set(struc_from_top.bond_types)
            assert set(struc.angle_types) == set(struc_from_top.angle_types)
            assert set(struc.dihedral_types) == set(struc_from_top.dihedral_types)
            assert set(struc.rb_torsion_types) == set(struc_from_top.rb_torsion_types)

            # Detail comparisions
            for i in range(len(struc.atoms)):
                assert struc_from_top.atoms[i].name == struc.atoms[i].name
                assert struc_from_top.atoms[i].atom_type == struc.atoms[i].atom_type

            for i in range(len(struc.bonds)):
                assert struc_from_top.bonds[i].atom1.name == struc.bonds[i].atom1.name
                assert struc_from_top.bonds[i].atom2.name == struc.bonds[i].atom2.name
                assert struc_from_top.bonds[i].type == struc.bonds[i].type

            for i in range(len(struc.angles)):
                assert struc_from_top.angles[i].atom1.name == struc.angles[i].atom1.name
                assert struc_from_top.angles[i].atom2.name == struc.angles[i].atom2.name
                assert struc_from_top.angles[i].atom3.name == struc.angles[i].atom3.name
                assert struc_from_top.angles[i].type == struc.angles[i].type

            for i in range(len(struc.dihedrals)):
                assert struc_from_top.dihedrals[i].atom1.name == struc.dihedrals[i].atom1.name
                assert struc_from_top.dihedrals[i].atom2.name == struc.dihedrals[i].atom2.name
                assert struc_from_top.dihedrals[i].atom3.name == struc.dihedrals[i].atom3.name
                assert struc_from_top.dihedrals[i].atom4.name == struc.dihedrals[i].atom4.name
                assert struc_from_top.dihedrals[i].type == struc.dihedrals[i].type

