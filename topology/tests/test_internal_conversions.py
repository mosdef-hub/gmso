import pytest

from topology.tests.base_test import BaseTest
from topology.core.site import Site
from topology.core.dihedral import Dihedral
from topology.core.dihedral_type import DihedralType
from topology.lib.potential_templates import RyckaertBellemansTorsionPotential
from topology.lib.potential_templates import OPLSTorsionPotential
from topology.utils.conversions import convert_ryckaert_to_opls
from topology.utils.conversions import convert_opls_to_ryckaert

import unyt as u
import numpy as np

class TestInternalConversions(BaseTest):
    def test_rb_to_opls(self):

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

        connection_type = DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=params)

        connection_members = [Site(),Site(),Site(),Site()]

        # Create connections
        rb_connection = Dihedral(
                connection_members=connection_members,
                connection_type=connection_type)

        # The connection type for this dihedral will be replaced
        # upon the call to convert_ryckaert_to_opls
        opls_connection = Dihedral(
                connection_members=connection_members,
                connection_type=connection_type)

        # Convert connection to OPLS
        convert_ryckaert_to_opls(opls_connection)

        # Pick some angles to check
        angles = [ -2.38, -1.31, -0.44, 0.0, 0.26, 0.92, 1.84, 3.10 ]

        for angle in angles:
            rb_val = eval_rb(rb_connection,angle)
            opls_val = eval_opls(opls_connection,angle)
            assert np.isclose(rb_val.value,
                              opls_val.value)


    def test_opls_to_rb(self):

        params = { 'k0' : 1.38   * u.Unit('kJ/mol'),
                   'k1' : -0.51  * u.Unit('kJ/mol'),
                   'k2' : 2.2    * u.Unit('kJ/mol'),
                   'k3' : -0.25  * u.Unit('kJ/mol'),
                   'k4' : 1.44   * u.Unit('kJ/mol')
                 }

        name = OPLSTorsionPotential().name
        expression = OPLSTorsionPotential().expression
        variables = OPLSTorsionPotential().independent_variables

        connection_type = DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=params)

        connection_members = [Site(),Site(),Site(),Site()]

        # Create connections
        opls_connection = Dihedral(
                connection_members=connection_members,
                connection_type=connection_type)


        # The connection type for this dihedral will be replaced
        # upon the call to convert_opls_to_ryckaert
        rb_connection = Dihedral(
                connection_members=connection_members,
                connection_type=connection_type)

        # Convert connection to OPLS
        convert_opls_to_ryckaert(rb_connection)

        # Pick some angles to check
        angles = [ -2.38, -1.31, -0.44, 0.0, 0.26, 0.92, 1.84, 3.10 ]

        for angle in angles:
            rb_val = eval_rb(rb_connection,angle)
            opls_val = eval_opls(opls_connection,angle)
            assert np.isclose(rb_val.value,
                              opls_val.value)


    def eval_rb(connection,angle):
        """Evaluate the RB torsion at a given angle (radians)"""

        params = connection.connection_type.parameters
        rb = 0.0
        for n in range(6):
            rb += params['c'+str(n)]*np.cos(angle-np.pi)**n
        return rb

    def eval_opls(connection,angle):
        """Evaluate the OPLS torsion at a given angle (radians)"""

        params = connection.connection_type.parameters
        opls = (1./2. * params['k0'] +
                1./2. * params['k1'] * (1. + np.cos(angle)) +
                1./2. * params['k2'] * (1. - np.cos(2.*angle)) +
                1./2. * params['k3'] * (1. + np.cos(3.*angle)) +
                1./2. * params['k4'] * (1. - np.cos(4.*angle)))
        return opls


