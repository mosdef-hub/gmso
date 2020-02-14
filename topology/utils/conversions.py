
import topology as topo
from topology.lib.potential_templates import RyckaertBellemansTorsionPotential
from topology.exceptions import TopologyError

def convert_opls_to_ryckaert(connection):
    """Convert an OPLS dihedral to Ryckaert-Bellemans dihedral

    Equations taken from:
        http://manual.gromacs.org/documentation/2019/
        reference-manual/functions/bonded-interactions.html
    """

    if connection.connection_type.parameters['k0'] != 0:
        raise TopologyError('Cannot convert from OPLS to Ryckaert'
                '-Bellemans dihedrals if "k0" != 0')

    f1 = connection.connection_type.parameters['k1']
    f2 = connection.connection_type.parameters['k2']
    f3 = connection.connection_type.parameters['k3']
    f4 = connection.connection_type.parameters['k4']

    new_params = {
            'c0' : (f2 + 0.5 * (f1 + f3)),
            'c1' : (0.5 * (-f1 + 3 * f3)),
            'c2' : (-f2 + 4 * f4),
            'c3' : (-2 * f3),
            'c4' : (-4 * f4),
            'c5' : 0.0 * u.Unit('kJ/mol')
    }

    name = RyckaertBellemansTorsionPotential().name
    expression = RyckaertBellemansTorsionPotential().expression
    variables = RyckaertBellemansTorsionPotential().independent_variables

    new_connection_type = topo.DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=new_params)

    connection.connection_type = new_connection_type


