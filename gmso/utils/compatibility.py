import sympy

from gmso.exceptions import EngineIncompatibilityError


def check_compatibility(topology, accepted_potentials):
    """Compare the potentials in a topology against a list of accepted potential templates"""
    for atom_type in topology.atom_types:
        if not _check_single_potential(atom_type, accepted_potentials):
            raise EngineIncompatibilityError

    for connection_type in topology.connection_types:
        if not _check_single_potential(connection_type, accepted_potentials):
            raise EngineIncompatibilityError


def _check_single_potential(potential, accepted_potentials):
    """Checks to see if a single given potential is in the list of accepted potentials"""
    for ref in accepted_potentials:
        if ref.independent_variables == potential.independent_variables:
            if sympy.simplify(ref.expression - potential.expression) == 0:
                return True
    return False
