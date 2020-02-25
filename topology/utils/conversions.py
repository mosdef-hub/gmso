
import sympy
import unyt as u

import topology as topo
from topology.lib.potential_templates import RyckaertBellemansTorsionPotential
from topology.lib.potential_templates import OPLSTorsionPotential
from topology.exceptions import TopologyError

def convert_opls_to_ryckaert(opls_connection_type):
    """Convert an OPLS dihedral to Ryckaert-Bellemans dihedral

    Equations taken/modified from:
        http://manual.gromacs.org/documentation/2019/
        reference-manual/functions/bonded-interactions.html

    NOTE: the conventions defining the dihedral angle are different
    for OPLS and RB torsions. OPLS torsions are defined with
    phi_cis = 0 while RB torsions are defined as phi_trans = 0.
    """

    valid_connection_type = False
    if ( opls_connection_type.independent_variables ==
         OPLSTorsionPotential().independent_variables ):
        if sympy.simplify(opls_connection_type.expression -
                          OPLSTorsionPotential().expression) == 0:
            valid_connection_type = True
    if not valid_connection_type:
        raise TopologyError('Cannot use convert_opls_to_ryckaert '
            'function to convert a ConnectionType that is not an '
            'OPLSTorsionPotential')

    f0 = opls_connection_type.parameters['k0']
    f1 = opls_connection_type.parameters['k1']
    f2 = opls_connection_type.parameters['k2']
    f3 = opls_connection_type.parameters['k3']
    f4 = opls_connection_type.parameters['k4']

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

    ryckaert_connection_type = topo.DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=converted_params)

    return ryckaert_connection_type

def convert_ryckaert_to_opls(ryckaert_connection_type):
    """Convert Ryckaert-Bellemans dihedral to OPLS

    NOTE: the conventions defining the dihedral angle are different
    for OPLS and RB torsions. OPLS torsions are defined with
    phi_cis = 0 while RB torsions are defined as phi_trans = 0.
    """

    valid_connection_type = False
    if ( ryckaert_connection_type.independent_variables ==
         RyckaertBellemansTorsionPotential().independent_variables ):
        if sympy.simplify(ryckaert_connection_type.expression -
                RyckaertBellemansTorsionPotential().expression) == 0:
            valid_connection_type = True
    if not valid_connection_type:
        raise TopologyError('Cannot use convert_ryckaert_to_opls '
            'function to convert a ConnectionType that is not an '
            'RyckaertBellemansTorsionPotential')


    c0 = ryckaert_connection_type.parameters['c0']
    c1 = ryckaert_connection_type.parameters['c1']
    c2 = ryckaert_connection_type.parameters['c2']
    c3 = ryckaert_connection_type.parameters['c3']
    c4 = ryckaert_connection_type.parameters['c4']
    c5 = ryckaert_connection_type.parameters['c5']

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

    opls_connection_type = topo.DihedralType(
            name=name,
            expression=expression,
            independent_variables=variables,
            parameters=converted_params)

    return opls_connection_type

