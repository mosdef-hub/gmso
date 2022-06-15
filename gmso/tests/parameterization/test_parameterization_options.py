import random

import forcefield_utilities as ffutils
import pytest

from gmso.core.forcefield import ForceField
from gmso.parameterization.parameterize import apply
from gmso.parameterization.topology_parameterizer import (
    GMSOParameterizationError,
)
from gmso.tests.parameterization.parameterization_base_test import (
    ParameterizationBaseTest,
)


class TestParameterizationOptions(ParameterizationBaseTest):
    def test_parameterization_error_different_scaling_factors(
        self, ethane_methane_top
    ):
        ff1 = ForceField()
        ff1.name = "FF1"
        ff1.scaling_factors = {
            "electrostatic14Scale": 1.0,
            "columbic14Scale": 2.0,
        }
        ff2 = ForceField()
        ff2.name = "FF2"
        ff2.scaling_factors = {
            "electrostatic14Scale": 3.0,
            "columbic14Scale": 2.0,
        }

        with pytest.raises(GMSOParameterizationError):
            apply(ethane_methane_top, {"Ethane": ff1, "Methane": ff2})

    def test_parameterization_different_combining_rule(
        self, ethane_methane_top
    ):
        ff1 = ForceField()
        ff1.name = "FF1"
        ff1.scaling_factors = {
            "electrostatic14Scale": 1.0,
            "columbic14Scale": 1.0,
        }
        ff1.combining_rule = "lorrentz"
        ff2 = ForceField()
        ff2.name = "FF2"
        ff2.scaling_factors = {
            "electrostatic14Scale": 1.0,
            "columbic14Scale": 1.0,
        }

        ff2.combining_rule = "geometric"

        with pytest.raises(GMSOParameterizationError):
            apply(ethane_methane_top, {"Ethane": ff1, "Methane": ff2})

    def test_different_ffs_apply(self, ethane_methane_top):
        opls = ffutils.FoyerFFs().load(ffname="oplsaa").to_gmso_ff()
        ethane_methane_top.identify_connections()
        apply(ethane_methane_top, {"Ethane": opls, "Methane": opls})
        assert ethane_methane_top.combining_rule == "geometric"
        for key, v in opls.scaling_factors.items():
            assert ethane_methane_top.scaling_factors[key] == v

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
