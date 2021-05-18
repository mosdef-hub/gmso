import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

from gmso.external.convert_parmed import from_parmed
from gmso.formats.gro import read_gro, write_gro
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn, has_parmed, import_

if has_parmed:
    pmd = import_("parmed")


@pytest.mark.skipif(not has_parmed, reason="ParmEd is not installed")
class TestGro(BaseTest):
    def test_read_gro(self):
        top = read_gro(get_fn("acn.gro"))

        assert top.name == "ACN"
        assert top.n_sites == 6
        assert_allclose_units(
            top.box.lengths, 4 * np.ones(3) * u.nm, rtol=1e-5, atol=1e-8
        )

        top = read_gro(get_fn("350-waters.gro"))

        assert top.name == "Generic title"
        assert top.n_sites == 1050
        assert_allclose_units(
            top.box.lengths, 2.20866 * np.ones(3) * u.nm, rtol=1e-5, atol=1e-8
        )

    def test_wrong_n_atoms(self):
        with pytest.raises(ValueError):
            read_gro(get_fn("too_few_atoms.gro"))
        with pytest.raises(ValueError):
            read_gro(get_fn("too_many_atoms.gro"))

    def test_write_gro(self):
        top = from_parmed(pmd.load_file(get_fn("ethane.gro"), structure=True))

        write_gro(top, "out.gro")

    def test_write_gro_non_orthogonal(self):
        top = from_parmed(pmd.load_file(get_fn("ethane.gro"), structure=True))
        top.box.angles = u.degree * [90, 90, 120]

        write_gro(top, "out.gro")
