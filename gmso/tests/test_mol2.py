import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

from gmso.external.convert_mol2 import from_mol2
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn


class TestMol2(BaseTest):
    def test_read_mol2(self):
        top = from_mol2(get_fn("parmed.mol2"))

        assert top.name == "parmed"
        assert top.n_sites == 8
        assert_allclose_units(
            top.box.lengths,
            ([8.2693, 7.9100, 6.6460] * u.Å).to("nm"),
            rtol=1e-5,
            atol=1e-8,
        )

        top = from_mol2(get_fn("tip3p.mol2"))

        assert top.name == "tip3p"
        assert top.n_sites == 3
        assert_allclose_units(
            top.box.lengths, 3.0130 * np.ones(3) * u.Å, rtol=1e-5, atol=1e-8
        )

        top = from_mol2(get_fn("vmd.mol2"))

        assert top.name == "vmd"
        assert top.n_sites == 6
        assert len(top.bonds) == 5
        assert top.bonds[0].connection_members[0] == top.sites[0]
        assert top.box == None

    def test_wrong_path(self):
        with pytest.raises(OSError):
            from_mol2(get_fn("not_a_file.mol2"))
        top = from_mol2(get_fn("ethane.gro"))
        assert len(top.sites) == 0
        assert len(top.bonds) == 0

    def test_broken_files(self):
        match = "The record type indicator @<TRIPOS>MOLECULE_extra_text\n is not supported"
        with pytest.warns(UserWarning) as record:
            from_mol2(get_fn("broken.mol2"))
        assert record[0].message.args[0] == match
        # with pytest.raises(UserWarning):
        #    from_mol2(get_fn("broken.mol2"))
