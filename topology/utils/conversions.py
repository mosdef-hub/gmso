
import topology as topo
from topology.lib.potential_templates import RyckaertBellemansTorsionPotential
from topology.lib.potential_templates import OPLSTorsionPotential
from topology.exceptions import TopologyError

import unyt as u

def convert_opls_to_ryckaert(connection):
    """Convert an OPLS dihedral to Ryckaert-Bellemans dihedral

    Equations taken/modified from:
        http://manual.gromacs.org/documentation/2019/
        reference-manual/functions/bonded-interactions.html
    """

    f0 = connection.connection_type.parameters['k0']
    f1 = connection.connection_type.parameters['k1']
    f2 = connection.connection_type.parameters['k2']
    f3 = connection.connection_type.parameters['k3']
    f4 = connection.connection_type.parameters['k4']

    converted_params = {
            'c0' : (f2 + 0.5 * (f0 + f1 + f3)),
            'c1' : (0.5 * (-f1 + 3. * f3)),
            'c2' : (-f2 + 4. * f4),
            'c3' : (-2. * f3),
            'c4' : (-4. * f4),
            'c5' : 0. * u.Unit('kJ/mol')
    }

    name = RyckaertBellemansTorsionPotential().name
    expression = RyckaertBellemansTorsionPotential().expression
    variables = RyckaertBellemansTorsionPotential().independent_variables

    updated_connection_type = topo.DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=converted_params)

    connection.connection_type = updated_connection_type

def convert_ryckaert_to_opls(connection):
    """Convert Ryckaert-Bellemans dihedrals to OPLS"""

    c0 = connection.connection_type.parameters['c0']
    c1 = connection.connection_type.parameters['c1']
    c2 = connection.connection_type.parameters['c2']
    c3 = connection.connection_type.parameters['c3']
    c4 = connection.connection_type.parameters['c4']
    c5 = connection.connection_type.parameters['c5']

    if c5 != 0.0:
        raise TopologyError('Cannot convert Ryckaert-Bellemans dihedral '
                'to OPLS dihedral if c5 is not equal to zero.')

    converted_params = {
            'k0' : 2. * (c0 + c1 + c2 + c3 + c4),
            'k1' : (-2. * c1 - (3./2.) * c3),
            'k2' : (-c2 - c4),
            'k3' : ((-1./2.) * c3),
            'k4' : ((-1./4.) * c4)
    }

    name = OPLSTorsionPotential().name
    expression = OPLSTorsionPotential().expression
    variables = OPLSTorsionPotential().independent_variables

    updated_connection_type = topo.DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=converted_params)

    connection.connection_type = updated_connection_type

