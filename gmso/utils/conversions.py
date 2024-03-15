"""Module for standard conversions needed in molecular simulations."""

import re
from functools import lru_cache

import numpy as np
import symengine
import sympy
import unyt as u
from unyt.dimensions import length, mass, time

import gmso
from gmso.exceptions import GMSOError
from gmso.lib.potential_templates import PotentialTemplateLibrary

templates = PotentialTemplateLibrary()


@lru_cache(maxsize=128)
def _constant_multiplier(pot1, pot2):
    # TODO: Doc string
    # TODO: Test outputs
    # TODO: Check speed
    try:
        constant = symengine.expand(pot1.expression / pot2.expression)
        # constant = sympy.simplify(pot1.expression / pot2.expression)
        if constant.is_number:
            for eq_term in pot1.expression.args:
                if eq_term.is_symbol:
                    key = str(eq_term)
                    return {key: pot1.parameters[key] * float(constant)}
    except Exception:
        # return nothing if the sympy conversion errors out
        pass
    return None


sympy_conversionsList = [_constant_multiplier]


@lru_cache(maxsize=128)
def _try_sympy_conversions(pot1, pot2):
    # TODO: Doc string
    # TODO: Test outputs
    # TODO: Check speed
    convertersList = []
    for conversion in sympy_conversionsList:
        convertersList.append(conversion(pot1, pot2))
    completed_conversions = np.where(convertersList)[0]
    if len(completed_conversions) > 0:  # check to see if any conversions worked
        return convertersList[
            completed_conversions[0]
        ]  # return first completed value
    return None


def convert_topology_expressions(top, expressionMap={}):
    """Convert from one parameter form to another.

    Parameters
    ----------
    expressionMap : dict, default={}
        map with keys of the potential type and the potential to change to

    Examples
    --------
    Convert from RB torsions to OPLS torsions
    top.convert_expressions({"dihedrals": "OPLSTorsionPotential"})
    """
    # TODO: Raise errors

    # Apply from predefined conversions or easy sympy conversions
    conversions_map = {
        (
            "OPLSTorsionPotential",
            "RyckaertBellemansTorsionPotential",
        ): convert_opls_to_ryckaert,
        (
            "RyckaertBellemansTorsionPotential",
            "OPLSTorsionPotential",
        ): convert_ryckaert_to_opls,
        (
            "RyckaertBellemansTorsionPotential",
            "FourierTorsionPotential",
        ): convert_ryckaert_to_opls,
    }  # map of all accessible conversions currently supported

    for conv in expressionMap:
        # check all connections with these types for compatibility
        for conn in getattr(top, conv):
            current_expression = getattr(conn, conv[:-1] + "_type")
            if (
                current_expression.name == expressionMap[conv]
            ):  # check to see if we can skip this one
                # TODO: Do something instead of just comparing the names
                continue

            # convert it using pre-defined conversion functions
            conversion_from_conversion_toTuple = (
                current_expression.name,
                expressionMap[conv],
            )
            if (
                conversion_from_conversion_toTuple in conversions_map
            ):  # Try mapped conversions
                new_conn_type = conversions_map.get(
                    conversion_from_conversion_toTuple
                )(current_expression)
                setattr(conn, conv[:-1] + "_type", new_conn_type)
                continue

            # convert it using sympy expression conversion
            new_potential = templates[expressionMap[conv]]
            modified_connection_parametersDict = _try_sympy_conversions(
                current_expression, new_potential
            )
            if modified_connection_parametersDict:  # try sympy conversions
                current_expression.name = new_potential.name
                current_expression.expression = new_potential.expression
                current_expression.parameters.update(
                    modified_connection_parametersDict
                )

    return top


@lru_cache(maxsize=128)
def convert_opls_to_ryckaert(opls_connection_type):
    """Convert an OPLS dihedral to Ryckaert-Bellemans dihedral.

    Equations taken/modified from:
        http://manual.gromacs.org/documentation/2019/
        reference-manual/functions/bonded-interactions.html

    NOTE: the conventions defining the dihedral angle are different
    for OPLS and RB torsions. OPLS torsions are defined with
    phi_cis = 0 while RB torsions are defined as phi_trans = 0.
    """
    # TODO: this function really converts the fourier torsion to rb, not opls
    opls_torsion_potential = templates["FourierTorsionPotential"]
    valid_connection_type = False
    if (
        opls_connection_type.independent_variables
        == opls_torsion_potential.independent_variables
    ):
        if (
            sympy.simplify(
                opls_connection_type.expression
                - opls_torsion_potential.expression
            )
            == 0
        ):
            valid_connection_type = True
    if not valid_connection_type:
        raise GMSOError(
            "Cannot use convert_opls_to_ryckaert "
            "function to convert a ConnectionType that is not an "
            "OPLSTorsionPotential"
        )

    f0 = opls_connection_type.parameters["k0"]
    f1 = opls_connection_type.parameters["k1"]
    f2 = opls_connection_type.parameters["k2"]
    f3 = opls_connection_type.parameters["k3"]
    f4 = opls_connection_type.parameters["k4"]

    converted_params = {
        "c0": (f2 + 0.5 * (f0 + f1 + f3)),
        "c1": (0.5 * (-f1 + 3.0 * f3)),
        "c2": (-f2 + 4.0 * f4),
        "c3": (-2.0 * f3),
        "c4": (-4.0 * f4),
        "c5": 0.0 * u.Unit("kJ/mol"),
    }
    ryckaert_bellemans_torsion_potential = templates[
        "RyckaertBellemansTorsionPotential"
    ]
    name = ryckaert_bellemans_torsion_potential.name
    expression = ryckaert_bellemans_torsion_potential.expression
    variables = ryckaert_bellemans_torsion_potential.independent_variables

    ryckaert_connection_type = gmso.DihedralType(
        name=name,
        expression=expression,
        independent_variables=variables,
        parameters=converted_params,
        member_types=opls_connection_type.member_types,
    )

    return ryckaert_connection_type


@lru_cache(maxsize=128)
def convert_ryckaert_to_opls(ryckaert_connection_type):
    """Convert Ryckaert-Bellemans dihedral to OPLS.

    NOTE: the conventions defining the dihedral angle are different
    for OPLS and RB torsions. OPLS torsions are defined with
    phi_cis = 0 while RB torsions are defined as phi_trans = 0.
    """
    fourier_connection_type = convert_ryckaert_to_fourier(
        ryckaert_connection_type
    )
    opls_torsion_potential = templates["OPLSTorsionPotential"]
    converted_params = {
        k: fourier_connection_type.parameters.get(k, None)
        for k in ["k1", "k2", "k3", "k4"]
    }

    name = opls_torsion_potential.name
    expression = opls_torsion_potential.expression
    variables = opls_torsion_potential.independent_variables

    opls_connection_type = gmso.DihedralType(
        name=name,
        expression=expression,
        independent_variables=variables,
        parameters=converted_params,
        member_types=ryckaert_connection_type.member_types,
    )

    return opls_connection_type


@lru_cache(maxsize=128)
def convert_ryckaert_to_fourier(ryckaert_connection_type):
    """Convert Ryckaert-Bellemans dihedral to Fourier.

    NOTE: the conventions defining the dihedral angle are different
    for OPLS and RB torsions. OPLS torsions are defined with
    phi_cis = 0 while RB torsions are defined as phi_trans = 0.
    """
    templates = PotentialTemplateLibrary()
    ryckaert_bellemans_torsion_potential = templates[
        "RyckaertBellemansTorsionPotential"
    ]
    fourier_torsion_potential = templates["FourierTorsionPotential"]

    valid_connection_type = False
    if (
        ryckaert_connection_type.independent_variables
        == ryckaert_bellemans_torsion_potential.independent_variables
    ):
        if (
            sympy.simplify(
                ryckaert_connection_type.expression
                - ryckaert_bellemans_torsion_potential.expression
            )
            == 0
        ):
            valid_connection_type = True
    if not valid_connection_type:
        raise GMSOError(
            "Cannot use convert_ryckaert_to_fourier "
            "function to convert a ConnectionType that is not an "
            "RyckaertBellemansTorsionPotential"
        )

    c0 = ryckaert_connection_type.parameters["c0"]
    c1 = ryckaert_connection_type.parameters["c1"]
    c2 = ryckaert_connection_type.parameters["c2"]
    c3 = ryckaert_connection_type.parameters["c3"]
    c4 = ryckaert_connection_type.parameters["c4"]
    c5 = ryckaert_connection_type.parameters["c5"]

    if c5 != 0.0:
        raise GMSOError(
            "Cannot convert Ryckaert-Bellemans dihedral "
            "to Fourier dihedral if c5 is not equal to zero."
        )

    converted_params = {
        "k0": 2.0 * (c0 + c1 + c2 + c3 + c4),
        "k1": (-2.0 * c1 - (3.0 / 2.0) * c3),
        "k2": (-c2 - c4),
        "k3": ((-1.0 / 2.0) * c3),
        "k4": ((-1.0 / 4.0) * c4),
    }

    name = fourier_torsion_potential.name
    expression = fourier_torsion_potential.expression
    variables = fourier_torsion_potential.independent_variables

    fourier_connection_type = gmso.DihedralType(
        name=name,
        expression=expression,
        independent_variables=variables,
        parameters=converted_params,
        member_types=ryckaert_connection_type.member_types,
    )

    return fourier_connection_type


def convert_kelvin_to_energy_units(
    energy_input_unyt,
    energy_output_unyt_units_str,
):
    """Convert the Kelvin (K) energy unit to a standard energy unit.

    Check to see if the unyt energy value is in Kelvin (K) and converts it to
    another energy unit (Ex: kcal/mol, kJ/mol, etc.).  Otherwise, it passes thru the
    existing unyt values.

    Parameters
    ----------
    energy_input_unyt : unyt.unyt_quantity
        The energy in units of 'energy / other_units' Example: 'energy/mol/angstroms**2' or 'K/angstroms**2'.
        NOTE: The only valid temperature unit for thermal energy is Kelvin (K).
    energy_output_unyt_units_str : str (unyt valid energy units only)
        (Example - 'kcal/mol', 'kJ/mol', or any '(length)**2*(mass)/(time)**2' , but not 'K')
        The energy units which a Kelvin (K) energy is converted into.
        It does not convert any energy unit if the the energy_input_unyt is not in Kelvin (K).
        NOTE and WARNING: the energy units of kcal, kJ will be accepted due to the way the
        Unyt module does not accept mol as a recorded unit; however, this will result in
        incorrect results from the intended function.


    Returns
    -------
    energy_output_unyt : unyt.unyt_quantity
        If the energy_input_unyt is in Kelvin (K), it converted to the specified energy_output_unyt_units_str.
        Otherwise, it passes through the existing unyt values.

    """
    # check for input errors
    # if not isinstance(energy_input_unyt, type(u.unyt_quantity(1, "K"))):
    if not isinstance(energy_input_unyt, u.unyt_quantity):
        print_error_message = (
            f"ERROR: The entered energy_input_unyt value is a {type(energy_input_unyt)}, "
            f"not a {type(u.Kelvin)}."
        )
        raise ValueError(print_error_message)

    if not isinstance(energy_output_unyt_units_str, str):
        print_error_message = (
            f"ERROR: The entered energy_output_unyt_units_str value is a {type(energy_output_unyt_units_str)}, "
            f"not a {type('string')}."
        )
        raise ValueError(print_error_message)

    # check for K energy units and convert them to normal energy units;
    # otherwise, just pass thru the original unyt units
    if energy_output_unyt_units_str in ["K"]:
        print_error_message = f"ERROR: The entered energy_output_unyt_units_str can not be in K energy units."
        raise ValueError(print_error_message)

    elif (length) ** 2 * (mass) / (time) ** 2 != u.unyt_quantity(
        1, energy_output_unyt_units_str
    ).units.dimensions:
        print_error_message = (
            f"ERROR: The entered energy_output_unyt_units_str value must be in units of energy/mol, "
            f"(length)**2*(mass)/(time)**2, but not in K energy units."
        )
        raise ValueError(print_error_message)

    if ("K") in str(energy_input_unyt.units) and "temperature" in str(
        energy_input_unyt.units.dimensions
    ):
        K_to_energy_conversion_constant = u.unyt_quantity(1, "K").to_value(
            energy_output_unyt_units_str, equivalence="thermal"
        )
        energy_output_unyt = (
            energy_input_unyt
            / u.Kelvin
            * u.unyt_quantity(
                K_to_energy_conversion_constant, energy_output_unyt_units_str
            )
        )
    else:
        energy_output_unyt = energy_input_unyt

    return energy_output_unyt


def convert_params_units(
    potentials, expected_units_dim, base_units, ref_values
):
    """Convert parameters' units in the potential to that specified in the base_units."""
    converted_potentials = list()
    for potential in potentials:
        converted_params = dict()
        for parameter in potential.parameters:
            unit_dim = expected_units_dim[parameter]
            ind_units = re.sub("[^a-zA-Z]+", " ", unit_dim).split()
            for unit in ind_units:
                if unit != "angle":
                    unit_dim = unit_dim.replace(unit, f"{base_units[unit]}")
                else:
                    # angle doesn't show up, degree or radian does
                    unit_dim = unit_dim.replace(unit, str(base_units[unit]))

            converted_params[parameter] = potential.parameters[parameter].to(
                u.Unit(unit_dim, registry=base_units.registry)
            )
        potential.parameters = converted_params
        converted_potentials.append(potential)
    return converted_potentials
