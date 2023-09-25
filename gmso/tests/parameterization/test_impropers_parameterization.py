from pathlib import Path

import pytest
import unyt as u

from gmso.core.topology import Topology
from gmso.core.views import PotentialFilters
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.parameterization.parameterize import apply
from gmso.parameterization.topology_parameterizer import ParameterizationError
from gmso.tests.parameterization.parameterization_base_test import (
    ParameterizationBaseTest,
)
from gmso.tests.utils import get_path


class TestImpropersParameterization(ParameterizationBaseTest):
    def test_improper_parameterization(self, fake_improper_ff_gmso, ethane):
        ethane.identify_connections()
        apply(ethane, fake_improper_ff_gmso, ignore_params=list())

        lib = PotentialTemplateLibrary()
        template_improper_type = lib["PeriodicImproperPotential"]

        assert (
            len(
                ethane.improper_types(
                    filter_by=PotentialFilters.UNIQUE_NAME_CLASS
                )
            )
            == 2
        )
        for improper_type in ethane.improper_types:
            assert improper_type.expression == template_improper_type.expression
            assert improper_type.member_classes in {
                ("CT", "CT", "HC", "HC"),
                ("CT", "HC", "HC", "HC"),
            }

            if improper_type.member_classes == ("CT", "CT", "HC", "HC"):
                assert u.allclose_units(
                    improper_type.parameters["phi_eq"], [180.0] * u.degree
                )
                assert u.allclose_units(
                    improper_type.parameters["k"], [4.6024] * u.kJ / u.mol
                )
                assert u.allclose_units(
                    improper_type.parameters["n"], [2] * u.dimensionless
                )

            elif improper_type.member_classes == ("CT", "CT", "HC", "HC"):
                assert u.allclose_units(
                    improper_type.parameters["phi_eq"], [180.0] * u.degree
                )
                assert u.allclose_units(
                    improper_type.parameters["k"], [2.5560] * u.kJ / u.mol
                )
                assert u.allclose_units(
                    improper_type.parameters["n"], [2] * u.dimensionless
                )

    def test_improper_assertion_error(self, ethane_methane_top, oplsaa_gmso):
        with pytest.raises(ParameterizationError):
            apply(ethane_methane_top, oplsaa_gmso, ignore_params=list())

    @pytest.mark.parametrize(
        "mol2_loc",
        [
            get_path("methyl_benzene_aa.mol2"),
            get_path("benzene_aa.mol2"),
            get_path("ethyl_benzene_aa.mol2"),
        ],
        ids=lambda p: Path(p).stem,
    )
    def test_benzene_aa_ff(self, mol2_loc, benzene_alkane_aa_ff_gmso):
        gmso_top = Topology.load(filename=mol2_loc)
        apply(gmso_top, benzene_alkane_aa_ff_gmso, identify_connections=True)
        for improper in gmso_top.impropers:
            if improper.improper_type:
                if [
                    improper.member_types[0],
                    set(improper.member_types[1:]),
                ] == ["CH_sp2", {"HCE", "CH_sp2", "C_sp2"}]:
                    params = improper.improper_type.get_parameters()
                    assert u.allclose_units(params["k"], 4.0 * u.kJ / u.mol)
                    assert u.allclose_units(params["n"], 2.0 * u.dimensionless)
                    assert u.allclose_units(params["phi_eq"], 0 * u.radian)

                elif [
                    improper.member_types[0],
                    set(improper.member_types[1:]),
                ] == ["C_sp2", {"CH3_sp3", "CH_sp2", "CH_sp2"}]:
                    params = improper.improper_type.get_parameters()
                    assert u.allclose_units(params["k"], 4.0 * u.kJ / u.mol)
                    assert u.allclose_units(params["n"], 1.0 * u.dimensionless)
                    assert u.allclose_units(
                        params["phi_eq"], 180 * u.degree, atol=10e-4
                    )
            else:
                assert len(improper.member_classes) == 4
                assert set(improper.member_classes) not in [
                    {"CE", "HCE", "CE", "CE"},
                    {"CE", "CT", "CE", "CE"},
                ]
