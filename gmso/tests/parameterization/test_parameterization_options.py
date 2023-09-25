import random
from copy import deepcopy

import forcefield_utilities as ffutils
import mbuild as mb
import pytest
from foyer.exceptions import FoyerError

import gmso
from gmso.core.forcefield import ForceField
from gmso.core.topology import Topology
from gmso.external.convert_mbuild import from_mbuild
from gmso.parameterization.parameterize import apply
from gmso.parameterization.topology_parameterizer import ParameterizationError
from gmso.tests.parameterization.parameterization_base_test import (
    ParameterizationBaseTest,
)
from gmso.tests.utils import get_path


class TestParameterizationOptions(ParameterizationBaseTest):
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

    def test_populating_member_types(self, ethane):
        ethane.identify_connections()
        opls = ffutils.FoyerFFs().load(ffname="oplsaa").to_gmso_ff()
        apply(top=ethane, forcefields=opls, remove_untyped=True)
        for connection in ethane.connections:
            connection_type = connection.connection_type
            assert (
                connection_type.member_types and connection_type.member_classes
            )
            for i in range(len(connection_type.member_classes)):
                assert (
                    opls.atom_types[connection_type.member_types[i]].atomclass
                    == connection_type.member_classes[i]
                )

    def test_different_ffs_apply(self, ethane_methane_top):
        opls = ffutils.FoyerFFs().load(ffname="oplsaa").to_gmso_ff()
        opls_copy = ffutils.FoyerFFs().load(ffname="oplsaa").to_gmso_ff()
        opls_copy.scaling_factors = {
            "nonBonded14Scale": 1.2,
            "electrostatics14Scale": 1.5,
        }
        ethane_methane_top.identify_connections()
        apply(
            ethane_methane_top,
            {"Ethane": opls, "Methane": opls_copy},
            "molecule",
        )

        assert ethane_methane_top.combining_rule == "geometric"
        assert (
            ethane_methane_top.get_lj_scale(
                molecule_id="Ethane", interaction="14"
            )
            == opls.scaling_factors["nonBonded14Scale"]
            == 0.5
        )
        assert (
            ethane_methane_top.get_electrostatics_scale(
                molecule_id="Ethane", interaction="14"
            )
            == opls.scaling_factors["electrostatics14Scale"]
            == 0.5
        )

        assert (
            ethane_methane_top.get_lj_scale(
                molecule_id="Methane", interaction="14"
            )
            == opls_copy.scaling_factors["nonBonded14Scale"]
            == 1.2
        )
        assert (
            ethane_methane_top.get_electrostatics_scale(
                molecule_id="Methane", interaction="14"
            )
            == opls_copy.scaling_factors["electrostatics14Scale"]
            == 1.5
        )

    """ Change to molecule"""

    def test_no_molecule_dict_ff(self, oplsaa_gmso):
        top = Topology(name="topWithNoMolecule")
        with pytest.warns(UserWarning):
            apply(top, {"moleculeA": oplsaa_gmso})
        assert not top.is_typed()

    def test_missing_group_name_ff(self, oplsaa_gmso):
        top = Topology(name="top1")
        for j in range(0, 10, 2):
            top.add_site(gmso.Atom(name=f"Atom_{j+1}", group="groupB"))
        with pytest.warns(
            UserWarning,
            match=r"Group/molecule groupB will not be parameterized, as the forcefield "
            r"to parameterize it is missing.",
        ):
            apply(top, {"groupA": oplsaa_gmso}, match_ff_by="group")

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

    @pytest.mark.parametrize(
        "speedup_by_molgraph, speedup_by_moltag",
        [(False, False), (True, False), (False, True), (True, True)],
    )
    def test_speedup_options(
        self,
        ethane_box_with_methane,
        oplsaa_gmso,
        speedup_by_molgraph,
        speedup_by_moltag,
    ):
        ethane_box_with_methane.identify_connections()
        apply(
            ethane_box_with_methane,
            oplsaa_gmso,
            identify_connections=False,
            speedup_by_molgraph=speedup_by_molgraph,
            speedup_by_moltag=speedup_by_moltag,
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

    def test_remove_untyped(self, oplsaa_gmso):
        isopropane = mb.load("C(C)C", smiles=True)
        top1 = gmso.external.from_mbuild(isopropane)
        top1.identify_connections()
        assert top1.n_impropers != 0
        apply(top1, oplsaa_gmso, remove_untyped=False)
        assert top1.n_impropers != 0

        top2 = gmso.external.from_mbuild(isopropane)
        top2.identify_connections()
        assert top2.n_impropers != 0
        apply(top2, oplsaa_gmso, remove_untyped=True)
        assert top2.n_impropers == 0

    def test_match_ff_by_molecule(self, ethane_box_with_methane, oplsaa_gmso):
        ethane_box_with_methane.identify_connections()
        ff_dict = {"Ethane": oplsaa_gmso, "Methane": oplsaa_gmso}
        apply(
            ethane_box_with_methane,
            ff_dict,
            match_ff_by="molecule",
            identify_connections=False,
            speedup_by_molgraph=True,
            speedup_by_moltag=True,
        )
        assert ethane_box_with_methane.atom_types is not None

    def test_match_ff_by_group(self, ethane_box_with_methane, oplsaa_gmso):
        ethane_box_with_methane.identify_connections()
        for site in ethane_box_with_methane.sites:
            site.group = "Alkane"
        ff_dict = {
            "Alkane": oplsaa_gmso,
        }
        apply(
            ethane_box_with_methane,
            ff_dict,
            match_ff_by="group",
            identify_connections=False,
            speedup_by_molgraph=True,
            speedup_by_moltag=True,
        )
        assert ethane_box_with_methane.atom_types is not None

    @pytest.mark.parametrize(
        "speedup_by_molgraph, speedup_by_moltag, match_ff_by",
        [
            (False, False, "group"),
            (True, False, "group"),
            (False, True, "group"),
            (True, True, "group"),
            (False, False, "molecule"),
            (True, False, "molecule"),
            (False, True, "molecule"),
            (True, True, "molecule"),
        ],
    )
    def test_hierarchical_mol_structure(
        self,
        oplsaa_gmso,
        hierarchical_top,
        speedup_by_molgraph,
        speedup_by_moltag,
        match_ff_by,
    ):
        top = deepcopy(hierarchical_top)
        # Load forcefield dicts
        spce = ForceField(get_path("spce.xml"))
        if match_ff_by == "molecule":
            ff_dict = {
                "polymer": oplsaa_gmso,
                "cyclopentane": oplsaa_gmso,
                "water": spce,
            }
        elif match_ff_by == "group":
            ff_dict = {"sol1": oplsaa_gmso, "sol2": spce}
        else:
            raise ValueError("Unexpected value provided match_ff_by.")

        apply(
            top,
            ff_dict,
            speedup_by_molgraph=speedup_by_molgraph,
            speedup_by_moltag=speedup_by_moltag,
            match_ff_by=match_ff_by,
        )
