from gmso.formats.itp import read_itp
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path


class Testitp(BaseTest):
    def test_itp_LIQ(self):
        top = read_itp(get_path("LIQ.itp"))
        assert top is not None
        assert len(top.atom_types) == 14
        assert len(set([atype.name for atype in top.atom_types])) == 4
        assert len(top.bond_types) == 13
        # assert top.bonds[0].bond_type.parameters["k"] == 0.1529
        empty_set = set()
        for bond in top.bonds:
            b_type_tuple = (
                float(bond.bond_type.parameters["k"].value),
                float(bond.bond_type.parameters["r_eq"].value),
            )
            empty_set.add(b_type_tuple)
        assert len(empty_set) == 4

        for angle in top.angles:
            a_type_tuple = (
                float(angle.angle_type.parameters["k"].value),
                float(angle.angle_type.parameters["theta_eq"].value),
            )
            empty_set.add(a_type_tuple)
        assert len(empty_set) == 10

        for dihedral in top.dihedrals:
            d_type_tuple = (
                float(dihedral.dihedral_type.parameters["k"].value),
                float(dihedral.dihedral_type.parameters["phi_eq"].value),
                float(dihedral.dihedral_type.parameters["n"].value),
            )
            empty_set.add(d_type_tuple)
        assert len(empty_set) == 14

    def test_itp_PNB(self):
        top = read_itp(get_path("PNB.itp"))
        assert top is not None

        assert len(top.atom_types) == 1416
        assert len(set([atype.name for atype in top.atom_types])) == 9

        assert len(top.bond_types) == 1445
        # assert top.bonds[0].bond_type.parameters["k"] == 0.1529
        empty_set = set()
        for bond in top.bonds:
            b_type_tuple = (
                float(bond.bond_type.parameters["k"].value),
                float(bond.bond_type.parameters["r_eq"].value),
            )
            empty_set.add(b_type_tuple)
        assert len(empty_set) == 6

        for angle in top.angles:
            a_type_tuple = (
                float(angle.angle_type.parameters["k"].value),
                float(angle.angle_type.parameters["theta_eq"].value),
            )
            empty_set.add(a_type_tuple)
        assert len(empty_set) == 16

        for dihedral in top.dihedrals:
            d_type_tuple = (
                float(dihedral.dihedral_type.parameters["c0"].value),
                float(dihedral.dihedral_type.parameters["c1"].value),
                float(dihedral.dihedral_type.parameters["c2"].value),
                float(dihedral.dihedral_type.parameters["c3"].value),
                float(dihedral.dihedral_type.parameters["c4"].value),
                float(dihedral.dihedral_type.parameters["c5"].value),
            )
            empty_set.add(d_type_tuple)
        assert len(empty_set) == 28
