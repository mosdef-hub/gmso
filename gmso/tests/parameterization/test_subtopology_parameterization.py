import forcefield_utilities as ffutils
import mbuild as mb
import pytest
from mbuild.lib.molecules import Ethane, Methane

from gmso.external.convert_mbuild import from_mbuild
from gmso.parameterization.parameterize import apply
from gmso.tests.base_test import BaseTest


class TestSubTopologyParameterization(BaseTest):
    @pytest.fixture(scope="session")
    def ethane_methane_top(self):
        cmpd = mb.Compound()
        cmpd.add(Ethane())
        cmpd.add(Methane())
        gmso_top = from_mbuild(cmpd)
        gmso_top.identify_connections()
        return gmso_top

    def test_different_ffs_apply(self, ethane_methane_top):
        opls = ffutils.FoyerFFs().load(ffname="oplsaa").to_gmso_ff()
        apply(ethane_methane_top, {"Ethane": opls, "Methane": opls})
