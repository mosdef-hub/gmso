import random

import forcefield_utilities as ffutils
import mbuild as mb
import pytest
from mbuild.lib.molecules import Ethane, Methane

from gmso.external.convert_mbuild import from_mbuild
from gmso.parameterization.parameterize import apply
from gmso.tests.base_test import BaseTest


class TestParameterizationOptions(BaseTest):
    @pytest.fixture(scope="session")
    def ethane_methane_top(self):
        cmpd = mb.Compound()
        cmpd.add(Ethane())
        cmpd.add(Methane())
        gmso_top = from_mbuild(cmpd)
        gmso_top.identify_connections()
        return gmso_top

    @pytest.fixture(scope="session")
    def ethane_box_with_methane(self):
        cmpd_box = mb.fill_box([Ethane(), Methane()], [50, 50], density=1.0)
        return from_mbuild(cmpd_box)

    def test_different_ffs_apply(self, ethane_methane_top):
        opls = ffutils.FoyerFFs().load(ffname="oplsaa").to_gmso_ff()
        ethane_methane_top.identify_connections()
        apply(ethane_methane_top, {"Ethane": opls, "Methane": opls})

    def test_isomporhic_speedups(self, ethane_box_with_methane, oplsaa_gmso):
        ethane_box_with_methane.identify_connections()
        apply(
            ethane_box_with_methane,
            oplsaa_gmso,
            identify_connections=False,
            identify_connected_components=True,
        )

        ethane_subtops = list(
            filter(
                lambda subtop: subtop.name == "Ethane",
                ethane_box_with_methane.subtops,
            )
        )
        methane_subtops = list(
            filter(
                lambda subtop: subtop.name == "Methane",
                ethane_box_with_methane.subtops,
            )
        )
        ethane_a = random.choice(ethane_subtops)
        ethane_b = random.choice(ethane_subtops)
        for atom_a, atom_b in zip(ethane_a.sites, ethane_b.sites):
            assert atom_a.atom_type == atom_b.atom_type
            assert atom_a.atom_type is not None

        methane_a = random.choice(methane_subtops)
        methane_b = random.choice(methane_subtops)
        for atom_a, atom_b in zip(methane_a.sites, methane_b.sites):
            assert atom_a.atom_type == atom_b.atom_type
            assert atom_a.atom_type is not None
