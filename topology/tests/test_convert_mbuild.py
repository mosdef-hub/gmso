import pytest

import unyt as u

from topology.external.convert_mbuild import from_mbuild, to_mbuild
from topology.tests.base_test import BaseTest
from topology.utils.io import has_mbuild
from topology.utils.testing import allclose


if has_mbuild:
    from mbuild.examples import Ethane
    from mbuild.box import Box

@pytest.mark.skipif(not has_mbuild, reason="mBuild is not installed")
class TestConvertMBuild(BaseTest):
    def test_from_mbuild(self):
        ethane = Ethane()
        top = from_mbuild(ethane)

        assert top.n_sites == 8
        assert top.n_connections == 7

    def test_full_conversion(self):
        ethane = Ethane()
        top = from_mbuild(ethane)

        new = to_mbuild(top)

        assert new.n_particles == 8
        assert new.n_bonds == 7

    def test_pass_box(self):
        ethane = Ethane()
        mb_box = Box(lengths=[3,3,3])

        top = from_mbuild(ethane, box=mb_box)
        assert allclose(top.box.lengths, [3,3,3]*u.nm)


    def test_pass_failed_box(self):
        ethane = Ethane()

        with pytest.raises(ValueError):
            top = from_mbuild(ethane, box=[3,3,3])

    def test_pass_box_periodicity(self):
        ethane = Ethane()
        ethane.periodicity = [2,2,2]
        top = from_mbuild(ethane)
        assert allclose(top.box.lengths, [2,2,2]*u.nm)

    def test_pass_box_bounding(self):
        ethane = Ethane()
        ethane.periodicity = [0,0,0]
        top = from_mbuild(ethane)
        assert allclose(top.box.lengths, 
                (ethane.boundingbox.lengths + [0.5, 0.5, 0.5]) * u.nm)

