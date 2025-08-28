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
        for bond in top.bonds:
            b_type_tuple = (
                float(bond.bond_type.parameters["k"].value),
                float(bond.bond_type.parameters["r_eq"].value),
            )
            empty_set.add(b_type_tuple)
            assert len(empty_set) == 4

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

        for bond in top.bonds:
            b_type_tuple = (
                float(bond.bond_type.parameters["k"].value),
                float(bond.bond_type.parameters["r_eq"].value),
            )
            empty_set.add(b_type_tuple)
            assert len(empty_set) == 6
