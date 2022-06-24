import random

import forcefield_utilities as ffutils
import pytest
from foyer.exceptions import FoyerError

import gmso
from gmso.core.forcefield import ForceField
from gmso.core.topology import Topology
from gmso.parameterization.parameterize import apply
from gmso.parameterization.topology_parameterizer import ParameterizationError
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

        with pytest.raises(ParameterizationError):
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

        with pytest.raises(ParameterizationError):
            apply(ethane_methane_top, {"Ethane": ff1, "Methane": ff2})

    def test_different_ffs_apply(self, ethane_methane_top):
        opls = ffutils.FoyerFFs().load(ffname="oplsaa").to_gmso_ff()
        ethane_methane_top.identify_connections()
        apply(ethane_methane_top, {"Ethane": opls, "Methane": opls})
        assert ethane_methane_top.combining_rule == "geometric"
        for key, v in opls.scaling_factors.items():
            assert ethane_methane_top.scaling_factors[key] == v

    """ Change to molecule"""

    def test_no_molecule_dict_ff(self, oplsaa_gmso):
        top = Topology(name="topWithNoMolecule")
        with pytest.raises(ParameterizationError):
            apply(top, {"moleculeA": oplsaa_gmso})

    def test_missing_group_name_ff(self, oplsaa_gmso):
        top = Topology(name="top1")
        for j in range(0, 10, 2):
            top.add_site(gmso.Atom(name=f"Atom_{j+1}", group="groupB"))
        with pytest.warns(
            UserWarning,
            match=r"Group groupB will not be parameterized, as the forcefield "
            r"to parameterize it is missing.",
        ):
            apply(top, {"groupA": oplsaa_gmso})

    def test_diff_combining_rules_error(self, ethane_methane_top):
        ff1 = ForceField()
        ff1.combining_rule = "lorentz"
        ff2 = ForceField()
        ff2.combining_rule = "geometric"
        with pytest.raises(ParameterizationError, match=""):
            apply(ethane_methane_top, {"Ethane": ff1, "Methane": ff2})

    def test_empty_ff_foyer_error(self, ethane_methane_top):
        with pytest.raises(FoyerError):
            apply(ethane_methane_top, ForceField())

    def test_empty_top_parameterization(self, oplsaa_gmso):
        with pytest.raises(FoyerError):
            apply(top=Topology(), forcefields=oplsaa_gmso)

    def test_isomorphic_speedups(self, ethane_box_with_methane, oplsaa_gmso):
        ethane_box_with_methane.identify_connections()
        apply(
            ethane_box_with_methane,
            oplsaa_gmso,
            identify_connections=False,
            identify_connected_components=True,
        )

        molecule_labels = ethane_box_with_methane.unique_site_labels("molecule")
        ethane_molecules = [
            label for label in molecule_labels if label.name == "Ethane"
        ]
        methane_molecules = [
            label for label in molecule_labels if label.name == "Methane"
        ]

        ethane_a = tuple(
            ethane_box_with_methane.iter_sites(
                "molecule", random.choice(ethane_molecules)
            )
        )
        ethane_b = tuple(
            ethane_box_with_methane.iter_sites(
                "molecule", random.choice(ethane_molecules)
            )
        )
        for atom_a, atom_b in zip(ethane_a, ethane_b):
            assert atom_a.atom_type == atom_b.atom_type
            assert atom_a.atom_type is not None

        methane_a = tuple(
            ethane_box_with_methane.iter_sites(
                "molecule", random.choice(methane_molecules)
            )
        )
        methane_b = tuple(
            ethane_box_with_methane.iter_sites(
                "molecule", random.choice(methane_molecules)
            )
        )
        for atom_a, atom_b in zip(methane_a, methane_b):
            assert atom_a.atom_type == atom_b.atom_type
            assert atom_a.atom_type is not None
