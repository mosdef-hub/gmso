from gmso.formats.itp import read_itp
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path


class Testitp(BaseTest):
    def test_itp(self):
        # top = Topology.load(get_fn("acn.gro"))
        top = read_itp(get_path("PNB.itp"))
        assert top != None

        assert len(top.atom_types) == 516
        assert len(set([atype.name for atype in top.atom_types])) == 6

        assert len(top.bond_types) == 545
        assert top.bonds[0].bond_type.parameters["k"] == 0.1529
        empty_set = set()
        for bond in top.bonds:
            b_type_tuple = (
                float(bond.bond_type.parameters["k"].value),
                float(bond.bond_type.parameters["r_eq"].value),
            )
            empty_set.add(b_type_tuple)
        assert len(empty_set) == 5

        for bond in top.bonds:
            b_type_tuple = (
                float(bond.bond_type.parameters["k"].value),
                float(bond.bond_type.parameters["r_eq"].value),
            )
            empty_set.add(b_type_tuple)
            assert len(empty_set) == 5
