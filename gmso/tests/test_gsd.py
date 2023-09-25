import gsd.hoomd
import mbuild as mb
import pytest
import unyt as u

from gmso.external.convert_mbuild import from_mbuild
from gmso.external.convert_parmed import from_parmed
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn, has_gsd, has_parmed, import_

if has_parmed:
    pmd = import_("parmed")


@pytest.mark.skipif(not has_gsd, reason="gsd is not installed")
@pytest.mark.skipif(not has_parmed, reason="ParmEd is not installed")
class TestGsd(BaseTest):
    def test_write_gsd_untyped(self):
        comp = mb.load("CCCC", smiles=True)
        system = mb.fill_box(comp, n_compounds=3, density=100)
        top = from_mbuild(system)
        top.identify_connections()
        top.save("out.gsd")
        with gsd.hoomd.open("out.gsd") as traj:
            snap = traj[0]
            assert all([i in snap.particles.types for i in ["C", "H"]])
            assert all([i in snap.bonds.types for i in ["C-C", "C-H"]])
            assert all([i in snap.angles.types for i in ["C-C-C", "C-C-H"]])
            assert all(
                [i in snap.dihedrals.types for i in ["C-C-C-C", "C-C-C-H"]]
            )

    def test_write_gsd(self, hierarchical_compound):
        top = from_mbuild(hierarchical_compound)
        top.save("out.gsd")

    def test_write_gsd_pmd(self):
        top = from_parmed(
            pmd.load_file(get_fn("ethane.top"), xyz=get_fn("ethane.gro"))
        )
        top.save("out.gsd")

    def test_write_gsd_non_orthogonal(self):
        top = from_parmed(
            pmd.load_file(get_fn("ethane.top"), xyz=get_fn("ethane.gro"))
        )
        top.box.angles = u.degree * [90, 90, 120]
        top.save("out.gsd")
