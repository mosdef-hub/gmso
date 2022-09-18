import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

from gmso import Topology
from gmso.external.convert_parmed import from_parmed
from gmso.formats.gro import read_gro, write_gro
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path
from gmso.utils.io import get_fn, has_mbuild, has_parmed, import_

if has_parmed:
    pmd = import_("parmed")


@pytest.mark.skipif(not has_parmed, reason="ParmEd is not installed")
class TestGro(BaseTest):
    def test_read_gro(self):
        top = Topology.load(get_fn("acn.gro"))

        assert top.name == "ACN"
        assert top.n_sites == 6
        assert_allclose_units(
            top.box.lengths, 4 * np.ones(3) * u.nm, rtol=1e-5, atol=1e-8
        )

        top = Topology.load(get_fn("350-waters.gro"))

        assert top.name == "Generic title"
        assert top.n_sites == 1050
        assert_allclose_units(
            top.box.lengths, 2.20866 * np.ones(3) * u.nm, rtol=1e-5, atol=1e-8
        )

    def test_wrong_n_atoms(self):
        with pytest.raises(IndexError):
            Topology.load(get_fn("too_few_atoms.gro"))
        with pytest.raises(ValueError):
            Topology.load(get_fn("too_many_atoms.gro"))

    def test_write_gro(self):
        top = from_parmed(pmd.load_file(get_fn("ethane.gro"), structure=True))
        top.save("out.gro")

    def test_write_gro_non_orthogonal(self):
        top = from_parmed(pmd.load_file(get_fn("ethane.gro"), structure=True))
        top.box.angles = u.degree * [90, 90, 120]
        top.save("out.gro")

    @pytest.mark.skipif(not has_mbuild, reason="mBuild not installed.")
    def test_benzene_gro(self):
        import mbuild as mb
        from mbuild.packing import fill_box

        from gmso.external import from_mbuild

        benzene = mb.load(get_fn("benzene.mol2"))
        benzene.children[0].name = "Benzene"
        box_of_benzene = fill_box(compound=benzene, n_compounds=5, density=1)
        top = from_mbuild(box_of_benzene)
        top.save("benzene.gro")

        reread = Topology.load("benzene.gro")
        for site, ref_site in zip(reread.sites, top.sites):
            assert site.molecule.name == ref_site.molecule.name[:3]
            assert site.molecule.number == ref_site.molecule.number

    def test_gro_read_molecule(self):
        top = Topology.load(get_path("benzene.gro"))
        for site in top.sites:
            assert site.molecule
            assert site.molecule.name == "Ben"
        assert len(top.unique_site_labels("molecule")) == 5
