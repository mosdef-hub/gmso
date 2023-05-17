import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

from gmso.core.dihedral_type import DihedralType
from gmso.exceptions import GMSOError
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.tests.base_test import BaseTest
from gmso.utils.conversions import (
    convert_opls_to_ryckaert,
    convert_ryckaert_to_fourier,
)


class TestInternalConversions(BaseTest):
    @pytest.fixture
    def templates(self):
        return PotentialTemplateLibrary()

    def test_invalid_connection_type(self, templates):
        params = {
            "c0": 1.53 * u.Unit("kJ/mol"),
            "c1": 0.76 * u.Unit("kJ/mol"),
            "c2": -0.22 * u.Unit("kJ/mol"),
            "c3": 3.55 * u.Unit("kJ/mol"),
            "c4": 0.94 * u.Unit("kJ/mol"),
            "c5": 0.0 * u.Unit("kJ/mol"),
        }
        ryckaert_bellemans_torsion_potential = templates[
            "RyckaertBellemansTorsionPotential"
        ]

        name = ryckaert_bellemans_torsion_potential.name
        expression = (
            str(ryckaert_bellemans_torsion_potential.expression) + " + 3 * psi"
        )
        variables = ["phi", "psi"]

        ryckaert_connection_type = DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=params,
        )

        with pytest.raises(GMSOError, match="Cannot use"):
            convert_ryckaert_to_fourier(ryckaert_connection_type)

        expression = "c0+c1+c2+c3+c4+c5+phi"
        variables = ryckaert_bellemans_torsion_potential.independent_variables
        ryckaert_connection_type = DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=params,
        )

        with pytest.raises(GMSOError, match="Cannot use"):
            convert_ryckaert_to_fourier(ryckaert_connection_type)

        # Pick some OPLS parameters at random
        params = {
            "k0": 1.38 * u.Unit("kJ/mol"),
            "k1": -0.51 * u.Unit("kJ/mol"),
            "k2": 2.2 * u.Unit("kJ/mol"),
            "k3": -0.25 * u.Unit("kJ/mol"),
            "k4": 1.44 * u.Unit("kJ/mol"),
        }

        opls_torsion_potential = templates["OPLSTorsionPotential"]
        name = opls_torsion_potential.name
        expression = (
            str(opls_torsion_potential.expression)
            + " + 0.5 * k1 * (1 + cos(psi))"
        )
        variables = ["phi", "psi"]

        opls_connection_type = DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=params,
        )

        with pytest.raises(GMSOError, match=""):
            ryckaert_connection_type = convert_opls_to_ryckaert(
                opls_connection_type
            )

        variables = opls_torsion_potential.independent_variables
        expression = "k0+k1+k2+k3+k4+phi"
        opls_connection_type = DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=params,
        )

        with pytest.raises(GMSOError, match=""):
            ryckaert_connection_type = convert_opls_to_ryckaert(
                opls_connection_type
            )

    def test_ryckaert_to_fourier(self, templates):
        # Pick some RB parameters at random
        params = {
            "c0": 1.53 * u.Unit("kJ/mol"),
            "c1": 0.76 * u.Unit("kJ/mol"),
            "c2": -0.22 * u.Unit("kJ/mol"),
            "c3": 3.55 * u.Unit("kJ/mol"),
            "c4": 0.94 * u.Unit("kJ/mol"),
            "c5": 0.0 * u.Unit("kJ/mol"),
        }

        ryckaert_bellemans_torsion_potential = templates[
            "RyckaertBellemansTorsionPotential"
        ]

        name = ryckaert_bellemans_torsion_potential.name
        expression = ryckaert_bellemans_torsion_potential.expression
        variables = ryckaert_bellemans_torsion_potential.independent_variables

        ryckaert_connection_type = DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=params,
        )

        # Convert connection to Fourier
        opls_connection_type = convert_ryckaert_to_fourier(
            ryckaert_connection_type
        )

        # Pick some angles to check
        angles = [-2.38, -1.31, -0.44, 0.0, 0.26, 0.92, 1.84, 3.10]

        for angle in angles:
            assert np.isclose(
                float(
                    ryckaert_connection_type.expression.subs(
                        [
                            (param, val)
                            for param, val in {
                                **ryckaert_connection_type.parameters,
                                "phi": angle - np.pi,
                            }.items()
                        ]
                    )
                ),
                float(
                    opls_connection_type.expression.subs(
                        [
                            (param, val)
                            for param, val in {
                                **opls_connection_type.parameters,
                                "phi": angle,
                            }.items()
                        ]
                    )
                ),
            )

    def test_opls_to_ryckaert(self, templates):
        # Pick some OPLS parameters at random
        params = {
            "k0": 1.38 * u.Unit("kJ/mol"),
            "k1": -0.51 * u.Unit("kJ/mol"),
            "k2": 2.2 * u.Unit("kJ/mol"),
            "k3": -0.25 * u.Unit("kJ/mol"),
            "k4": 1.44 * u.Unit("kJ/mol"),
        }

        opls_torsion_potential = templates["FourierTorsionPotential"]
        name = opls_torsion_potential.name
        expression = opls_torsion_potential.expression
        variables = opls_torsion_potential.independent_variables

        opls_connection_type = DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=params,
        )

        # Convert connection to RB
        ryckaert_connection_type = convert_opls_to_ryckaert(
            opls_connection_type
        )

        # Pick some angles to check
        angles = [-2.38, -1.31, -0.44, 0.0, 0.26, 0.92, 1.84, 3.10]

        for angle in angles:
            assert np.isclose(
                float(
                    ryckaert_connection_type.expression.subs(
                        [
                            (param, val)
                            for param, val in {
                                **ryckaert_connection_type.parameters,
                                "phi": angle - np.pi,
                            }.items()
                        ]
                    )
                ),
                float(
                    opls_connection_type.expression.subs(
                        [
                            (param, val)
                            for param, val in {
                                **opls_connection_type.parameters,
                                "phi": angle,
                            }.items()
                        ]
                    )
                ),
            )

    def test_double_conversion(self, templates):
        # Pick some OPLS parameters at random
        params = {
            "k0": 1.38 * u.Unit("kJ/mol"),
            "k1": -0.51 * u.Unit("kJ/mol"),
            "k2": 2.2 * u.Unit("kJ/mol"),
            "k3": -0.25 * u.Unit("kJ/mol"),
            "k4": 1.44 * u.Unit("kJ/mol"),
        }

        opls_torsion_potential = templates["FourierTorsionPotential"]

        name = opls_torsion_potential.name
        expression = opls_torsion_potential.expression
        variables = opls_torsion_potential.independent_variables

        opls_connection_type = DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=params,
        )

        # Convert connection to RB
        ryckaert_connection_type = convert_opls_to_ryckaert(
            opls_connection_type
        )

        # Convert connection back to OPLS
        final_connection_type = convert_ryckaert_to_fourier(
            ryckaert_connection_type
        )

        assert_allclose_units(
            [*opls_connection_type.parameters.values()],
            [*final_connection_type.parameters.values()],
            rtol=1e-5,
            atol=1e-8,
        )
