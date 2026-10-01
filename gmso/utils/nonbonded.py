"""Helpers for the non-bonded type pairs shared by the engine writers."""

import logging
from itertools import combinations_with_replacement

import numpy as np

from gmso.utils.expression import NullPotentialExpression
from gmso.utils.sorting import sort_by_types

logger = logging.getLogger(__name__)


def mix_sigma_epsilon(atom_types, combining_rule):
    """Return the combined sigma and epsilon of two atom types, keeping units.

    Parameters
    ----------
    atom_types : list-like of gmso.AtomType
        The two atom types of the pair to combine.
    combining_rule : str
        Either lorentz or geometric, as stored on Topology.combining_rule.

    Returns
    -------
    tuple of unyt.array.unyt_quantity
        The combined sigma and epsilon, in the units of the first atom type.

    Raises
    ------
    ValueError
        If combining_rule is neither lorentz nor geometric.
    """
    # Combine the magnitudes, not the quantities: unyt resolves mol against
    # Avogadro's number in a product, so kJ/mol * kJ/mol gives kJ**2.
    sigma_unit = atom_types[0].parameters["sigma"].units
    epsilon_unit = atom_types[0].parameters["epsilon"].units
    sigmas = [
        atom_type.parameters["sigma"].in_units(sigma_unit).value
        for atom_type in atom_types
    ]
    epsilons = [
        atom_type.parameters["epsilon"].in_units(epsilon_unit).value
        for atom_type in atom_types
    ]
    if combining_rule == "lorentz":
        sigma = np.mean(sigmas)
    elif combining_rule == "geometric":
        sigma = np.sqrt(sigmas[0] * sigmas[1])
    else:
        raise ValueError(f"Invalid combining rule provided ({combining_rule})")
    epsilon = np.sqrt(epsilons[0] * epsilons[1])
    return sigma * sigma_unit, epsilon * epsilon_unit


def explicit_pair_types(top, atom_type_names):
    """Map the sorted member types of every applicable PairPotentialType to it.

    Pair types naming an atom type outside the topology are left out, since
    apply() attaches every pair type in the force field. Insertion order is
    sorted by member types, so iterating the result is deterministic.

    Parameters
    ----------
    top : gmso.Topology
        The topology whose pairpotential_types to collect.
    atom_type_names : set of str
        Names of the atom types the writer is emitting.

    Returns
    -------
    dict
        Sorted member type tuple to gmso.PairPotentialType.
    """
    pair_types = {}
    for pairpotential_type in sorted(top.pairpotential_types, key=sort_by_types):
        members = sort_by_types(pairpotential_type)
        if not atom_type_names.issuperset(members):
            logger.debug(
                f"Pair potential type {pairpotential_type} is not written, because "
                f"the atom types {set(members) - atom_type_names} are not in the "
                f"topology {top}."
            )
            continue
        pair_types[members] = pairpotential_type
    return pair_types


def uncovered_null_pairs(atom_types, covered_pairs):
    """Return the bare atom types and the type pairs of theirs nothing covers.

    A bare atom type has no parameters to combine, so every pair it takes part
    in needs a PairPotentialType of its own.

    Parameters
    ----------
    atom_types : list-like of gmso.AtomType
        The atom types the writer is emitting.
    covered_pairs : container of tuple of str
        Sorted member type tuples that a PairPotentialType covers.

    Returns
    -------
    tuple of (list of str, list of tuple)
        Names of the bare atom types, and the type pairs touching one of them
        that covered_pairs does not contain. Both empty when none are bare.
    """
    null_names = {
        atom_type.name
        for atom_type in atom_types
        if isinstance(atom_type.potential_expression, NullPotentialExpression)
    }
    if not null_names:
        return [], []
    uncovered = [
        pair
        for pair in combinations_with_replacement(
            sorted(atom_type.name for atom_type in atom_types), 2
        )
        if null_names.intersection(pair) and pair not in covered_pairs
    ]
    return sorted(null_names), uncovered
