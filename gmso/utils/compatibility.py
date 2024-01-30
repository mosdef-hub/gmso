"""Determine if the parametrized gmso.topology can be written to an engine."""

from functools import lru_cache

import symengine
import sympy

from gmso.core.views import PotentialFilters
from gmso.exceptions import EngineIncompatibilityError


def check_compatibility(
    topology,
    accepted_potentials,
    site_pfilter=PotentialFilters.UNIQUE_NAME_CLASS,
    conn_pfilter=PotentialFilters.UNIQUE_ID,
):
    """
    Compare the potentials in a topology against a list of accepted potential templates.

    Parameters
    ----------
    topology: gmso.Topology
        The topology whose potentials to check.
    accepted_potentials: list
        A list of gmso.Potential objects to check against
    site_pfilter: gmso.core.view.PotentialFilter name, defaults=PotentialFilters.UNIQUE_NAME_CLASS
        A given name from the set of potential filters, or a user defined function that
        operates on each atom_type and returns the attributes of the atom_type to be
        considered as distinctive. In other words,
        e.g. site_pfilter = lambda atype: atype.name
    conn_pfilter: gmso.core.view.PotentialFilter name, default=PotentialFilters.UNIQUE_ID
        A given name from the set of potential filters, or a user defined function that
        operates on each connection_type and returns the attributes of the connection_type to be
        considered as distinctive.
        e.g. site_pfilter = lambda conn_type: conn_type.member_types

    Notes
    -----
    Pre-made potential identifiers that work for both site_pfilter and conn_pfilter.
        potential_identifiers = {
            PotentialFilters.UNIQUE_NAME_CLASS: get_name_or_class,
            PotentialFilters.UNIQUE_SORTED_NAMES: sort_by_types,
            PotentialFilters.UNIQUE_EXPRESSION: lambda p: str(p.expression),
            PotentialFilters.UNIQUE_PARAMETERS: get_parameters,
            PotentialFilters.UNIQUE_ID: lambda p: id(p),
            PotentialFilters.REPEAT_DUPLICATES: lambda _: str(uuid.uuid4()),
        }

    Returns
    -------
    potential_forms_dict: dict
        A dictionary mapping the atom_types and connection_types to the names
        of their matching template potential. Provides a standardized naming
        scheme for the potentials in the topology.

    """
    potential_forms_dict = dict()
    for atom_type in topology.atom_types(filter_by=site_pfilter):
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

    for connection_type in topology.connection_types(filter_by=conn_pfilter):
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


def _check_single_potential(potential, accepted_potentials):
    """Check to see if a single given potential is in the list of accepted potentials."""
    ind_var = potential.independent_variables
    u_dims = [
        symplify_str_eqn(para.units.dimensions)
        for para in potential.parameters.values()
    ]
    for i in range(len(u_dims)):
        if "temperature" in u_dims[i]:
            energy_dimsStr = str(replace_temp_with_energy(u_dims[i]))
            energy_dimsStr = symplify_str_eqn(energy_dimsStr)
            u_dims[i] = energy_dimsStr
    u_dims = set(u_dims)
    for ref in accepted_potentials:
        ref_ind_var = ref.independent_variables
        ref_u_dims = set(
            map(symplify_str_eqn, ref.expected_parameters_dimensions.values())
        )
        if len(ind_var) == len(ref_ind_var) and u_dims == ref_u_dims:
            if str(ref.expression) == str(potential.expression):
                return {potential: ref.name}
            else:
                if (
                    symengine.expand(ref.expression - potential.expression)
                    # sympy.simplify(ref.expression - potential.expression)
                    == 0
                ):
                    return {potential: ref.name}
    return False


def replace_temp_with_energy(dimsStr):
    """Track energy dimensions instead of temperature dimensions."""
    dimsStr = dimsStr.replace("(temperature)", "(length)**2*(mass)/(time)**2")
    return sympy.Symbol(dimsStr)


def symplify_str_eqn(eqnStr):
    """Take a string and use sympy to aggregate the expression."""
    if not isinstance(eqnStr, str):
        eqnStr = str(eqnStr)
    sym_eqn = symengine.expand(eqnStr)
    replace_symbolsSet = sym_eqn.free_symbols
    outStr = str(sym_eqn)
    for symbol in replace_symbolsSet:
        symbolStr = str(symbol)
        outStr = outStr.replace(symbolStr, f"({symbolStr})")
    return outStr
