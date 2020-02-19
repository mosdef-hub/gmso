
import topology as topo
import unyt as u
import numpy as np

from topology.lib.potential_templates import RyckaertBellemansTorsionPotential
from topology.lib.potential_templates import OPLSTorsionPotential
from topology.utils.conversions import convert_ryckaert_to_opls
from topology.utils.conversions import convert_opls_to_ryckaert

def main():

    # Check RB -> OPLS conversion
    for i in range(500):

        # Create an RB dihedral with random c0->c5 on (-5,5)
        params = {
                   'c0' : (10.*np.random.rand() - 5.0) * u.Unit('kJ/mol'),
                   'c1' : (10.*np.random.rand() - 5.0) * u.Unit('kJ/mol'),
                   'c2' : (10.*np.random.rand() - 5.0) * u.Unit('kJ/mol'),
                   'c3' : (10.*np.random.rand() - 5.0) * u.Unit('kJ/mol'),
                   'c4' : (10.*np.random.rand() - 5.0) * u.Unit('kJ/mol'),
                   'c5' : 0.0 * u.Unit('kJ/mol')
                }

        name = RyckaertBellemansTorsionPotential().name
        expression = RyckaertBellemansTorsionPotential().expression
        variables = RyckaertBellemansTorsionPotential().independent_variables

        connection_type = topo.DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=params)

        connection_members = [topo.Site(),topo.Site(),topo.Site(),topo.Site()]

        # Create connections
        rb_connection = topo.Dihedral(
                connection_members=connection_members,
                connection_type=connection_type)

        opls_connection = topo.Dihedral(
                connection_members=connection_members,
                connection_type=connection_type)

        # Convert connection to OPLS
        convert_ryckaert_to_opls(opls_connection)

        # Generate random angles to check
        angles = 2. * np.pi * np.random.rand(20) - np.pi

        for angle in angles:
            rb_val = eval_rb(rb_connection,angle)
            opls_val = eval_opls(opls_connection,angle)
            #print('Angle: {}\n  {}  {}'.format(angle,rb_val,opls_val))
            assert np.isclose(rb_val.value,
                              opls_val.value)


    # Check OPLS -> RB conversion
    for i in range(500):

        # Create an OPLS dihedral with random k0->k4 on (-50,50)
        params = {
                   'k0' : (10.*np.random.rand() - 5.0) * u.Unit('kJ/mol'),
                   'k1' : (10.*np.random.rand() - 5.0) * u.Unit('kJ/mol'),
                   'k2' : (10.*np.random.rand() - 5.0) * u.Unit('kJ/mol'),
                   'k3' : (10.*np.random.rand() - 5.0) * u.Unit('kJ/mol'),
                   'k4' : (10.*np.random.rand() - 5.0) * u.Unit('kJ/mol'),
                }

        name = OPLSTorsionPotential().name
        expression = OPLSTorsionPotential().expression
        variables = OPLSTorsionPotential().independent_variables

        connection_type = topo.DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=params)

        connection_members = [topo.Site(),topo.Site(),topo.Site(),topo.Site()]

        # Create connections
        rb_connection = topo.Dihedral(
                connection_members=connection_members,
                connection_type=connection_type)

        opls_connection = topo.Dihedral(
                connection_members=connection_members,
                connection_type=connection_type)

        # Convert connection to OPLS
        convert_opls_to_ryckaert(rb_connection)

        # Generate random angles to check
        angles = 2. * np.pi * np.random.rand(20) - np.pi

        for angle in angles:
            rb_val = eval_rb(rb_connection,angle)
            opls_val = eval_opls(opls_connection,angle)
            #print('Angle: {}\n  {}  {}'.format(angle,rb_val,opls_val))
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


if __name__ == "__main__":
    main()

