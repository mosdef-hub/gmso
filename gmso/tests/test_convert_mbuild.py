import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

import gmso
from gmso.core.atom import Atom
from gmso.core.topology import Topology as Top
from gmso.exceptions import GMSOError
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
        return mb.lib.molecules.Ethane()

    def test_from_mbuild_ethane(self, mb_ethane):
        import mbuild as mb

        top = from_mbuild(mb_ethane)

        assert top.n_sites == 8
        assert top.n_connections == 7
        for i in range(top.n_sites):
            assert isinstance(top.sites[i].element, gmso.Element)
            assert top.sites[i].name == top.sites[i].element.symbol
            assert top.sites[i].residue.name == "CH3"
            assert top.sites[i].molecule.name == "Ethane"

        unlabeled_top = from_mbuild(mb_ethane, parse_label=False)
        assert unlabeled_top.n_sites == 8
        assert unlabeled_top.n_connections == 7
        for site in unlabeled_top.sites:
            assert site.name == site.element.symbol
            assert site.residue is None
            assert site.molecule is None

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
            assert site.molecule[0] == "Ar"

    def test_from_mbuild_single_particle(self):
        compound = mb.Compound()
        top = from_mbuild(compound, parse_label=False)

        assert top.n_sites == 1
        assert top.n_connections == 0
        assert top.sites[0].residue == top.sites[0].molecule == None

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

        top = from_mbuild(top_cmpnd, parse_label=True)

        assert top.n_sites == 1
        assert top.sites[0].molecule == ("bot", 0)
        assert top.sites[0].residue == ("bot", 0)

    def test_4_layer_compound(self):
        l0_cmpnd = mb.Compound(name="l0")
        l1_cmpnd = mb.Compound(name="l1")
        l2_cmpnd = mb.Compound(name="l2")
        particle = mb.Compound(name="particle")

        l0_cmpnd.add(l1_cmpnd)
        l1_cmpnd.add(l2_cmpnd)
        l2_cmpnd.add(particle)

        l0_cmpnd.periodicity = [True, True, True]

        top = from_mbuild(l0_cmpnd, parse_label=True)

        assert top.n_sites == 1
        assert top.sites[0].molecule == ("particle", 0)

    def test_uneven_hierarchy(self):
        top_cmpnd = mb.Compound(name="top")
        mid_cmpnd = mb.Compound(name="mid")
        particle1 = mb.Compound(name="particle1")
        particle2 = mb.Compound(name="particle2")

        top_cmpnd.add(mid_cmpnd)
        top_cmpnd.add(particle1)
        mid_cmpnd.add(particle2)

        top_cmpnd.periodicity = [True, True, True]

        top = from_mbuild(top_cmpnd, parse_label=True)

        assert top.n_sites == 2
        for site in top.sites:
            if site.name == "particle2":
                assert site.group == "mid"
                assert site.molecule == ("particle2", 0)
            elif site.name == "particle1":
                assert site.molecule == ("particle1", 0)

    def test_pass_box(self, mb_ethane):
        mb_box = Box(lengths=[3, 3, 3])

        top = from_mbuild(mb_ethane, box=mb_box, parse_label=True)
        assert_allclose_units(
            top.box.lengths, [3, 3, 3] * u.nm, rtol=1e-5, atol=1e-8
        )

    def test_pass_failed_box(self, mb_ethane):
        with pytest.raises(ValueError):
            top = from_mbuild(mb_ethane, box=[3, 3, 3], parse_label=True)

    def test_pass_box_bounding(self, mb_ethane):
        mb_ethane.periodicity = [False, False, False]
        top = from_mbuild(mb_ethane, parse_label=True)
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

    def test_hierarchical_structure(self, hierarchical_top):
        for label in ("polymer", "water", "cyclopentane"):
            assert label in hierarchical_top.unique_site_labels(
                "molecule", name_only=True
            )
        for label in ("sol1", "sol2"):
            assert label in hierarchical_top.unique_site_labels(
                "group", name_only=True
            )

    @pytest.mark.skipif(not has_mbuild, reason="mBuild is not installed")
    def test_group_2_level_compound(self):
        mb_cpd = mb.Compound(name="_CH4", mass=12)
        filled_box = mb.fill_box(mb_cpd, n_compounds=2, density=0.01)
        filled_box.name = "group1"
        top = from_mbuild(filled_box)
        for site in top.sites:
            assert site.group == filled_box.name

    @pytest.mark.skipif(not has_mbuild, reason="mBuild is not installed")
    def test_custom_groups_from_compound(self):
        mb_cpd1 = mb.Compound(name="_CH4")

        first_bead = mb.Compound(name="_CH3")
        middle_bead = mb.Compound(name="_CH2")
        last_bead = mb.Compound(name="_CH3")
        mb_cpd2 = mb.Compound(name="Alkane")
        [mb_cpd2.add(cpd) for cpd in [first_bead, middle_bead, last_bead]]
        mb_cpd2.add_bond((first_bead, middle_bead))
        mb_cpd2.add_bond((last_bead, middle_bead))

        mb_cpd3 = mb.load("O", smiles=True)
        mb_cpd3.name = "O"

        filled_box1 = mb.fill_box(
            [mb_cpd1, mb_cpd2], n_compounds=[2, 2], box=[1, 1, 1]
        )
        filled_box1.name = "box1"
        filled_box2 = mb.fill_box(mb_cpd3, n_compounds=2, box=[1, 1, 1])
        filled_box2.name = "box2"

        top_box = mb.Compound()
        top_box.add(filled_box1)
        top_box.add(filled_box2)
        top_box.name = "top"

        list_of_groups = [
            (["top"], [14]),  # top level of hierarchy
            (["box1", "box2"], [8, 6]),  # middle level of hierarchy
            (["_CH4", "_CH2", "_CH3", "O"], [2, 2, 4, 6]),  # particle level
            (
                ["box2", "Alkane", "_CH4"],
                [6, 6, 2],
            ),  # multiple different levels
        ]
        for groups, n_groups in list_of_groups:
            top = from_mbuild(top_box, custom_groups=groups)
            assert np.all([site.group in groups for site in top.sites])
            for n, gname in zip(n_groups, groups):
                assert (
                    len([True for site in top.sites if site.group == gname])
                    == n
                )

    @pytest.mark.skipif(not has_mbuild, reason="mBuild is not installed")
    def test_single_custom_group(self):
        mb_cpd1 = mb.Compound(name="_CH4")
        mb_cpd2 = mb.Compound(name="_CH3")
        filled_box = mb.fill_box(
            [mb_cpd1, mb_cpd2], n_compounds=[2, 2], box=[1, 1, 1]
        )
        filled_box.name = "box1"

        top = from_mbuild(filled_box, custom_groups=filled_box.name)
        assert (
            len([True for site in top.sites if site.group == filled_box.name])
            == filled_box.n_particles
        )

    @pytest.mark.skipif(not has_mbuild, reason="mBuild is not installed")
    def test_bad_custom_groups_from_compound(self):
        mb_cpd1 = mb.Compound(name="_CH4")
        mb_cpd2 = mb.Compound(name="_CH3")
        filled_box = mb.fill_box(
            [mb_cpd1, mb_cpd2], n_compounds=[2, 2], box=[1, 1, 1]
        )

        with pytest.warns(Warning):
            top = from_mbuild(
                filled_box, custom_groups=["_CH4", "_CH3", "_CH5"]
            )

        with pytest.raises(GMSOError):
            top = from_mbuild(filled_box, custom_groups=["_CH4"])

        with pytest.raises(TypeError):
            top = from_mbuild(filled_box, custom_groups=mb_cpd1)

        with pytest.raises(TypeError):
            top = from_mbuild(filled_box, custom_groups=[mb_cpd1])

    @pytest.mark.skipif(not has_mbuild, reason="mBuild is not installed")
    def test_nontop_level_compound(self, mb_ethane):
        cpd = mb.Compound(name="top")
        cpd.add(mb_ethane)
        with pytest.raises(AssertionError):
            from_mbuild(mb_ethane)
