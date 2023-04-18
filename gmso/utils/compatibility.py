"""Determine if the parametrized gmso.topology can be written to an engine."""
from functools import lru_cache

import symengine
import sympy

from gmso.core.views import PotentialFilters
from gmso.exceptions import EngineIncompatibilityError


def check_compatibility(topology, accepted_potentials):
    """
    Compare the potentials in a topology against a list of accepted potential templates.

    Parameters
    ----------
    topology: gmso.Topology
        The topology whose potentials to check.
    accepted_potentials: list
        A list of gmso.Potential objects to check against

    Returns
    -------
    potential_forms_dict: dict
        A dictionary mapping the atom_types and connection_types to the names
        of their matching template potential. Provides a standardized naming
        scheme for the potentials in the topology.

    """
    potential_forms_dict = dict()
    for atom_type in topology.atom_types(
        filter_by=PotentialFilters.UNIQUE_NAME_CLASS
    ):
        potential_form = _check_single_potential(
            atom_type,
            accepted_potentials,
        )
        if not potential_form:
            raise EngineIncompatibilityError(
                f"Potential {atom_type} is not in the list of accepted_potentials {accepted_potentials}"
            )
        else:
            potential_forms_dict.update(potential_form)

    for connection_type in topology.connection_types(
        # filter_by=PotentialFilters.UNIQUE_NAME_CLASS
    ):
        potential_form = _check_single_potential(
            connection_type,
            accepted_potentials,
        )
        if not potential_form:
            raise EngineIncompatibilityError(
                f"Potential {connection_type} is not in the list of accepted_potentials {accepted_potentials}"
            )
        else:
            potential_forms_dict.update(potential_form)

    return potential_forms_dict


# @lru_cache()
def _check_single_potential(potential, accepted_potentials):
    """Check to see if a single given potential is in the list of accepted potentials."""
    ind_var = potential.independent_variables
    u_dims = {para.units.dimensions for para in potential.parameters.values()}
    for ref in accepted_potentials:
        ref_ind_var = ref.independent_variables
        ref_u_dims = set(ref.expected_parameters_dimensions.values())
        if len(ind_var) == len(ref_ind_var) and u_dims == ref_u_dims:
            if str(ref.expression) == str(potential.expression):
                return {potential: ref.name}
            else:
                print(symengine.expand(ref.expression - potential.expression))
                if (
                    symengine.expand(ref.expression - potential.expression)
                    # sympy.simplify(ref.expression - potential.expression)
                    == 0
                ):
                    return {potential: ref.name}
    return False
