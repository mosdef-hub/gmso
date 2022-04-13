import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

import gmso
from gmso.core.atom import Atom
from gmso.core.subtopology import SubTopology as SubTop
from gmso.core.topology import Topology as Top
from gmso.external.convert_mbuild import from_mbuild, to_mbuild
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn, has_mbuild

if has_mbuild:
    import mbuild as mb
    from mbuild.box import Box


@pytest.mark.skipif(not has_mbuild, reason="mBuild is not installed")
class TestConvertMBuild(BaseTest):
    @pytest.fixture
    def mb_ethane(self):
        return mb.load(get_fn("ethane.mol2"))

    def test_from_mbuild_ethane(self, mb_ethane):
        import mbuild as mb

        top = from_mbuild(mb_ethane)

        assert top.n_sites == 8
        assert top.n_connections == 7
        for i in range(top.n_sites):
            assert isinstance(top.sites[i].element, gmso.Element)
            assert top.sites[i].name == top.sites[i].element.symbol

        for site in top.sites:
            assert site.molecule == site.group == "CH3-CH3"

    def test_from_mbuild_argon(self, ar_system):
        # ar_system is a 3x3x3nm box filled with 100 argon sites using
        # mBuild, and then converted to topology via from_mbuild.

        top = ar_system

        assert top.n_sites == 100
        assert top.n_connections == 0
        for i in range(top.n_sites):
            assert isinstance(top.sites[i].element, gmso.Element)
            assert top.sites[i].name == top.sites[i].element.symbol

        for site in top.sites:
            assert site.molecule == site.group == "ar"

    def test_from_mbuild_single_particle(self):
        compound = mb.Compound()
        top = from_mbuild(compound)

        assert top.n_sites == 1
        assert top.n_connections == 0
        assert top.sites[0].molecule == "DefaultMolecule"
        assert top.sites[0].group == "DefaultGroup"

    def test_to_mbuild_name_none(self):
        top = Top()
        top.add_site(Atom(position=[0.0, 0.0, 0.0]))
        top.name = None
        compound = to_mbuild(top)

        assert compound.name == "Compound"

    def test_full_conversion(self, ethane):
        top = ethane

        new = to_mbuild(top)

        assert new.n_particles == 8
        assert new.n_bonds == 7
        for i in range(new.n_particles):
            assert np.isclose(new[i].xyz, top.sites[i].position.value).all()

    def test_3_layer_compound(self):
        top_cmpnd = mb.Compound(name="top")
        mid_cmpnd = mb.Compound(name="mid")
        bot_cmpnd = mb.Compound(name="bot")

        top_cmpnd.add(mid_cmpnd)
        mid_cmpnd.add(bot_cmpnd)

        top_cmpnd.periodicity = [True, True, True]

        top = from_mbuild(top_cmpnd)

        assert top.n_sites == 1
        assert top.sites[0].group == "mid"
        assert top.sites[0].molecule == "bot"

    def test_3_layer_top(self):
        top_top = Top()
        mid_top = SubTop()
        site = Atom(position=[0.0, 0.0, 0.0])

        top_top.add_subtopology(mid_top)
        mid_top.add_site(site)

        compound = to_mbuild(top_top)

        assert len(compound.children) == 1
        assert compound.children[0].n_particles == 1
        assert compound.n_particles == 1

    def test_4_layer_compound(self):
        l0_cmpnd = mb.Compound(name="l0")
        l1_cmpnd = mb.Compound(name="l1")
        l2_cmpnd = mb.Compound(name="l2")
        particle = mb.Compound(name="particle")

        l0_cmpnd.add(l1_cmpnd)
        l1_cmpnd.add(l2_cmpnd)
        l2_cmpnd.add(particle)

        l0_cmpnd.periodicity = [True, True, True]

        top = from_mbuild(l0_cmpnd)

        assert top.n_sites == 1
        assert top.sites[0].group == "l1"
        assert top.sites[0].molecule == "particle"

    def test_uneven_hierarchy(self):
        top_cmpnd = mb.Compound(name="top")
        mid_cmpnd = mb.Compound(name="mid")
        particle1 = mb.Compound(name="particle1")
        particle2 = mb.Compound(name="particle2")

        top_cmpnd.add(mid_cmpnd)
        top_cmpnd.add(particle1)
        mid_cmpnd.add(particle2)

        top_cmpnd.periodicity = [True, True, True]

        top = from_mbuild(top_cmpnd)

        assert top.n_sites == 2
        for site in top.sites:
            if site.name == "particle2":
                assert site.group == "mid"
                assert site.molecule == "particle2"
            elif site.name == "particle1":
                assert site.group == site.molecule == "particle1"

    def test_pass_box(self, mb_ethane):
        mb_box = Box(lengths=[3, 3, 3])

        top = from_mbuild(mb_ethane, box=mb_box)
        assert_allclose_units(
            top.box.lengths, [3, 3, 3] * u.nm, rtol=1e-5, atol=1e-8
        )

    def test_pass_failed_box(self, mb_ethane):
        with pytest.raises(ValueError):
            top = from_mbuild(mb_ethane, box=[3, 3, 3])

    def test_pass_box_bounding(self, mb_ethane):
        mb_ethane.periodicity = [False, False, False]
        top = from_mbuild(mb_ethane)
        assert_allclose_units(
            top.box.lengths,
            (mb_ethane.get_boundingbox().lengths) * u.nm,
            rtol=1e-5,
            atol=1e-8,
        )

    def test_empty_compound_name(self):
        compound = mb.load("CCOC", smiles=True)
        top = from_mbuild(compound)
        assert top.name is not None
