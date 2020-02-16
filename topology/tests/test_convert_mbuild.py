import pytest
import unyt as u

import topology
from topology.core.topology import Topology as Top
from topology.core.subtopology import SubTopology as SubTop
from topology.core.site import Site
from topology.external.convert_mbuild import from_mbuild, to_mbuild
from topology.tests.base_test import BaseTest
from topology.utils.io import get_fn, has_mbuild
from topology.utils.testing import allclose


if has_mbuild:
    import mbuild as mb
    from mbuild.box import Box

@pytest.mark.skipif(not has_mbuild, reason="mBuild is not installed")
class TestConvertMBuild(BaseTest):
    @pytest.fixture
    def ethane(self):
        return mb.load(get_fn('ethane.mol2'))

    def test_from_mbuild_ethane(self, ethane):
        import mbuild as mb
        top = from_mbuild(ethane)

        assert top.n_sites == 8
        assert top.n_subtops == 1
        assert top.subtops[0].n_sites == 8
        assert top.n_connections == 7
        for i in range(top.n_sites):
            assert isinstance(top.sites[i].element, topology.Element)
            assert top.sites[i].name == top.sites[i].element.symbol

    def test_from_mbuild_argon(self, ar_system):
        # ar_system is a 3x3x3nm box filled with 100 argon sites using
        # mBuild, and then converted to topology via from_mbuild.

        top = ar_system

        assert top.n_sites == 100
        assert top.n_subtops == 0
        assert top.n_connections == 0
        for i in range(top.n_sites):
            assert isinstance(top.sites[i].element, topology.Element)
            assert top.sites[i].name == top.sites[i].element.symbol

    def test_from_mbuild_single_particle(self):
        compound = mb.Compound()
        top = from_mbuild(compound)

        assert top.n_sites == 1
        assert top.n_subtops == 0
        assert top.n_connections == 0

    def test_to_mbuild_name_none(self):
        top = Top()
        top.add_site(Site())
        top.name = None
        compound = to_mbuild(top)

        assert compound.name == 'Compound'

    def test_full_conversion(self, ethane):
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

        top_cmpnd.periodicity = [1, 1, 1]

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

    def test_4_layer_compound(self):
        l0_cmpnd = mb.Compound()
        l1_cmpnd = mb.Compound()
        l2_cmpnd = mb.Compound()
        particle = mb.Compound()

        l0_cmpnd.add(l1_cmpnd)
        l1_cmpnd.add(l2_cmpnd)
        l2_cmpnd.add(particle)

        l0_cmpnd.periodicity = [1, 1, 1]

        top = from_mbuild(l0_cmpnd)

        assert top.n_sites == 1
        assert top.n_subtops == 1
        assert top.subtops[0].n_sites == 1
        assert top.subtops[0].sites[0] == top.sites[0]

    def test_pass_box(self, ethane):
        mb_box = Box(lengths=[3,3,3])

        top = from_mbuild(ethane, box=mb_box)
        assert allclose(top.box.lengths, [3,3,3]*u.nm)


    def test_pass_failed_box(self, ethane):
        with pytest.raises(ValueError):
            top = from_mbuild(ethane, box=[3,3,3])

    def test_pass_box_periodicity(self, ethane):
        ethane.periodicity = [2,2,2]
        top = from_mbuild(ethane)
        assert allclose(top.box.lengths, [2,2,2]*u.nm)

    def test_pass_box_bounding(self, ethane):
        ethane.periodicity = [0,0,0]
        top = from_mbuild(ethane)
        assert allclose(top.box.lengths,
                (ethane.boundingbox.lengths + [0.5, 0.5, 0.5]) * u.nm)
