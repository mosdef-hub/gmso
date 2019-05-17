import pytest

from topology.external.convert_mbuild import from_mbuild, to_mbuild
from topology.tests.base_test import BaseTest
from topology.utils.io import has_mbuild


if has_mbuild:
    from mbuild.examples import Ethane

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
