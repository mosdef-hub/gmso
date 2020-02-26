import pytest
import unyt as u
import numpy as np

from gmso.tests.base_test import BaseTest
from gmso.core.dihedral_type import DihedralType
from gmso.lib.potential_templates import RyckaertBellemansTorsionPotential
from gmso.lib.potential_templates import OPLSTorsionPotential
from gmso.utils.conversions import convert_ryckaert_to_opls
from gmso.utils.conversions import convert_opls_to_ryckaert
from gmso.exceptions import GMSOError

class TestInternalConversions(BaseTest):

    def test_invalid_connection_type(self):
        params = { 'c0' : 1.53  * u.Unit('kJ/mol'),
                   'c1' : 0.76  * u.Unit('kJ/mol'),
                   'c2' : -0.22 * u.Unit('kJ/mol'),
                   'c3' : 3.55  * u.Unit('kJ/mol'),
                   'c4' : 0.94  * u.Unit('kJ/mol'),
                   'c5' : 0.0   * u.Unit('kJ/mol')
                 }

        name = RyckaertBellemansTorsionPotential().name
        expression = RyckaertBellemansTorsionPotential().expression
        variables = ['phi','psi']

        ryckaert_connection_type = DihedralType(
                name=name,
                expression=expression,
                independent_variables=variables,
                parameters=params)

        with pytest.raises(GMSOError,match='Cannot use'):
            opls_connection_type = convert_ryckaert_to_opls(
                    ryckaert_connection_type)

        expression = 'c0+c1+c2+c3+c4+c5+phi'
        variables = RyckaertBellemansTorsionPotential().independent_variables
        ryckaert_connection_type = DihedralType(
                name=name,
                expression=expression,
                independent_variables=variables,
                parameters=params)

        with pytest.raises(GMSOError,match='Cannot use'):
            opls_connection_type = convert_ryckaert_to_opls(
                    ryckaert_connection_type)

        # Pick some OPLS parameters at random
        params = { 'k0' : 1.38   * u.Unit('kJ/mol'),
                   'k1' : -0.51  * u.Unit('kJ/mol'),
                   'k2' : 2.2    * u.Unit('kJ/mol'),
                   'k3' : -0.25  * u.Unit('kJ/mol'),
                   'k4' : 1.44   * u.Unit('kJ/mol')
                 }

        name = OPLSTorsionPotential().name
        expression = OPLSTorsionPotential().expression
        variables = ['phi','psi']

        opls_connection_type = DihedralType(
                name=name,
                expression=expression,
                independent_variables=variables,
                parameters=params)

        with pytest.raises(GMSOError,match=''):
            ryckaert_connection_type = convert_opls_to_ryckaert(
                    opls_connection_type)

        variables = OPLSTorsionPotential().independent_variables
        expression = 'k0+k1+k2+k3+k4+phi'
        opls_connection_type = DihedralType(
                name=name,
                expression=expression,
                independent_variables=variables,
                parameters=params)

        with pytest.raises(GMSOError,match=''):
            ryckaert_connection_type = convert_opls_to_ryckaert(
                    opls_connection_type)


    def test_ryckaert_to_opls(self):

        # Pick some RB parameters at random
        params = { 'c0' : 1.53  * u.Unit('kJ/mol'),
                   'c1' : 0.76  * u.Unit('kJ/mol'),
                   'c2' : -0.22 * u.Unit('kJ/mol'),
                   'c3' : 3.55  * u.Unit('kJ/mol'),
                   'c4' : 0.94  * u.Unit('kJ/mol'),
                   'c5' : 0.0   * u.Unit('kJ/mol')
                 }

        name = RyckaertBellemansTorsionPotential().name
        expression = RyckaertBellemansTorsionPotential().expression
        variables = RyckaertBellemansTorsionPotential().independent_variables

        ryckaert_connection_type = DihedralType(
                name=name,
                expression=expression,
                independent_variables=variables,
                parameters=params)

        # Convert connection to OPLS
        opls_connection_type = convert_ryckaert_to_opls(
                ryckaert_connection_type)

        # Pick some angles to check
        angles = [ -2.38, -1.31, -0.44, 0.0, 0.26, 0.92, 1.84, 3.10 ]

        for angle in angles:
            assert np.isclose(
                    float(ryckaert_connection_type.expression.subs(
                        [(param, val) for param, val in
                            {**ryckaert_connection_type.parameters,
                                'phi':angle-np.pi}.items()])),
                    float(opls_connection_type.expression.subs(
                        [(param, val) for param, val in
                            {**opls_connection_type.parameters,
                                'phi':angle}.items()])),
            )

    def test_opls_to_ryckaert(self):

        # Pick some OPLS parameters at random
        params = { 'k0' : 1.38   * u.Unit('kJ/mol'),
                   'k1' : -0.51  * u.Unit('kJ/mol'),
                   'k2' : 2.2    * u.Unit('kJ/mol'),
                   'k3' : -0.25  * u.Unit('kJ/mol'),
                   'k4' : 1.44   * u.Unit('kJ/mol')
                 }

        name = OPLSTorsionPotential().name
        expression = OPLSTorsionPotential().expression
        variables = OPLSTorsionPotential().independent_variables

        opls_connection_type = DihedralType(
                name=name,
                expression=expression,
                independent_variables=variables,
                parameters=params)

        # Convert connection to RB
        ryckaert_connection_type = convert_opls_to_ryckaert(
                opls_connection_type)

        # Pick some angles to check
        angles = [ -2.38, -1.31, -0.44, 0.0, 0.26, 0.92, 1.84, 3.10 ]

        for angle in angles:
            assert np.isclose(
                    float(ryckaert_connection_type.expression.subs(
                        [(param, val) for param, val in
                            {**ryckaert_connection_type.parameters,
                                'phi':angle-np.pi}.items()])),
                    float(opls_connection_type.expression.subs(
                        [(param, val) for param, val in
                            {**opls_connection_type.parameters,
                                'phi':angle}.items()])),
            )

    def test_double_conversion(self):

        # Pick some OPLS parameters at random
        params = { 'k0' : 1.38   * u.Unit('kJ/mol'),
                   'k1' : -0.51  * u.Unit('kJ/mol'),
                   'k2' : 2.2    * u.Unit('kJ/mol'),
                   'k3' : -0.25  * u.Unit('kJ/mol'),
                   'k4' : 1.44   * u.Unit('kJ/mol')
                 }

        name = OPLSTorsionPotential().name
        expression = OPLSTorsionPotential().expression
        variables = OPLSTorsionPotential().independent_variables

        opls_connection_type = DihedralType(
                name=name,
                expression=expression,
                independent_variables=variables,
                parameters=params)

        # Convert connection to RB
        ryckaert_connection_type = convert_opls_to_ryckaert(
                opls_connection_type)

        # Convert connection back to OPLS
        final_connection_type = convert_ryckaert_to_opls(
                ryckaert_connection_type)

        assert np.allclose([*opls_connection_type.parameters.values()],
                           [*final_connection_type.parameters.values()])
