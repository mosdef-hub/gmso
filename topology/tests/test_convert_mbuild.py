import pytest

from topology.core.topology import Topology as Top
from topology.core.subtopology import SubTopology as SubTop
from topology.core.site import Site
from topology.external.convert_mbuild import from_mbuild, to_mbuild
from topology.tests.base_test import BaseTest
from topology.utils.io import has_mbuild


if has_mbuild:
    import mbuild as mb
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

    def test_3_layer_compound(self):
        top_cmpnd = mb.Compound()
        mid_cmpnd = mb.Compound()
        bot_cmpnd = mb.Compound()

        top_cmpnd.add(mid_cmpnd)
        mid_cmpnd.add(bot_cmpnd)

        top = from_mbuild(top_cmpnd)

        assert top.n_sites == 1
        assert top.n_subtops == 1
        assert top.subtops[0].n_sites == 1

    def test_3_layer_top(self):
        top_top = Top()
        mid_top = SubTop()
        site = Site()

        top_top.add_subtopology(mid_top)
        mid_top.add_site(site)

        compound = to_mbuild(top_top)

        assert len(compound.children) == 1
        assert compound.children[0].n_particles == 1
        assert compound.n_particles == 1
