"""GMSO equation or expression comparisons."""

import os

import sympy
import unyt as u

from gmso.utils.io import get_fn


# compare Lennard-Jones (LJ) non-bonded equations
def evaluate_nonbonded_lj_format_with_scaler(new_lj_form, base_lj_form):
    """Compare a new Lennard-Jones (LJ) form to a base LJ form (new LJ form / base LJ form).

    If the new form is the same as the base form, other than a scaling factor,
    it labels the form and provides its scalar from the base form.

    Parameters
    ----------
    new_lj_form : str
        The new Lennard-Jones (LJ) form that will be compared or divided by the base form.
    base_lj_form : str
        The base LJ form, which is the standard form.

    Returns
    -------
    list, [form_output, form_scalar]
        form_output : str
            'LJ', if the new_lj_form variable is a LJ non-bonded form.
            None, if the new_lj_form variable is not a LJ non-bonded form.
        form_scalar : float
            float, if the new_lj_form variable is a LJ non-bonded form.
            None, if the new_lj_form variable is not a LJ non-bonded form.
    """
    try:
        (
            eqn_ratio,
            epsilon,
            sigma,
            r,
            Rmin,
            two,
        ) = sympy.symbols("eqn_ratio epsilon sigma r Rmin two")
        values = sympy.nonlinsolve(
            [
                eqn_ratio
                - sympy.sympify(new_lj_form) / sympy.sympify(base_lj_form),
                Rmin - sigma * two ** (1 / 6),
                two - 2,
            ],
            [eqn_ratio, Rmin, two],
        )

        form_scalar = float(list(values)[0][0])
        form_output = "LJ"

    except:
        form_scalar = None
        form_output = None

    return [form_output, form_scalar]


# compare Mie non-bonded equations
def evaluate_nonbonded_mie_format_with_scaler(new_mie_form, base_mie_form):
    """Compare a new Mie form to a base Mie form (new Mie form / base Mie form).

    If the new form is the same as the base form, other than a scaling factor,
    it labels the form and provides its scalar from the base form.

    Parameters
    ----------
    new_mie_form : str
        The new Mie form that will be compared or divided by the base form.
    base_mie_form : str
        The base Mie form, which is the standard form.

    Returns
    -------
    list, [form_output, form_scalar]
        form_output : str
            'Mie', if the new_mie_form variable is a Mie non-bonded form.
            None, if the new_mie_form variable is not a Mie non-bonded form.
        form_scalar : float
            float, if the new_mie_form variable is a Mie non-bonded form.
            None, if the new_mie_form variable is not a Mie non-bonded form.
    """
    try:
        (
            eqn_ratio,
            epsilon,
            sigma,
            r,
            n,
        ) = sympy.symbols("eqn_ratio epsilon sigma r n")
        values = sympy.nonlinsolve(
            [
                eqn_ratio
                - sympy.sympify(new_mie_form) / sympy.sympify(base_mie_form),
            ],
            [eqn_ratio],
        )

        form_scalar = float(list(values)[0][0])
        form_output = "Mie"

    except:
        form_scalar = None
        form_output = None

    return [form_output, form_scalar]


# compare Exp6 non-bonded equations
def evaluate_nonbonded_exp6_format_with_scaler(new_exp6_form, base_exp6_form):
    """Compare a new Exp6 form to a base Exp6 form (new Mie form / base Mie form).

    If the new form is the same as the base form, other than a scaling factor,
    it labels the form and provides its scalar from the base form.

    Parameters
    ----------
    new_exp6_form : str
        The new Exp6 form that will be compared or divided by the base form.
    base_exp6_form : str
        The base Exp6 form, which is the standard form.

    Returns
    -------
    list, [form_output, form_scalar]
        form_output : str
            'Exp6', if the new_exp6_form variable is an Exp6 non-bonded form.
            None, if the new_exp6_form variable is not an Exp6 non-bonded form.
        form_scalar : float
            float, if the new_exp6_form variable is an Exp6 non-bonded form.
            None, if the new_exp6_form variable is not an Exp6 non-bonded form.
    """
    try:
        (
            eqn_ratio,
            epsilon,
            sigma,
            r,
            Rmin,
            alpha,
        ) = sympy.symbols("eqn_ratio epsilon sigma r Rmin alpha")
        values = sympy.nonlinsolve(
            [
                eqn_ratio
                - sympy.sympify(new_exp6_form) / sympy.sympify(base_exp6_form),
            ],
            [eqn_ratio],
        )

        form_scalar = float(list(values)[0][0])
        form_output = "Exp6"

    except:
        form_scalar = None
        form_output = None

    return [form_output, form_scalar]


def get_atom_type_expressions_and_scalars(atom_types_dict):
    """Get the field file expressions and its scalar from the base forms used from the GMSO XML file.

    The expression base forms are as follows.
    'LJ' = "4*epsilon * ((Rmin/r)**12 - 2*(Rmin/r)**6)"
    'Mie' = "(n/(n-m)) * (n/m)**(m/(n-m)) * epsilon * ((sigma/r)**n - (sigma/r)**m)",
    'Exp6' = "epsilon*alpha/(alpha-6) * (6/alpha*exp(alpha*(1-r/Rmin)) - (Rmin/r)**6)"

    Note that when n=12 and m=6 in the 'Mie' form, it is the exact same as the 'LJ' form.
    Note that the 'LJ' form above is equivalent to "epsilon * ((Rmin/r)**12 - 2*(Rmin/r)**6)"

    Parameters
    ----------
    atom_types_dict: dict
        A nested dict with the all the residues as the main keys.

        {'residue_name': {
        'expression': expression_form_for_all_atom_types,
        'atom_types': Topology.atom_types}}
        }

        Example with an LJ expression

        atom_types_dict = {
        'ETH': {'expression': epsilon*(-sigma**6/r**6 + sigma**12/r**12),
        'atom_types': (<AtomType opls_135,
        expression: epsilon*(-sigma**6/r**6 + sigma**12/r**12),
        id: 140214732710736,
        atomclass: CT>, <AtomType opls_140,
        expression: epsilon*(-sigma**6/r**6 + sigma**12/r**12),
        id: 140214730466576,
        atomclass: HC>)}
        }

    Returns
    -------
    atom_types_data_expression_data_dict : dict
        The dictionary to append with the residue name and the GMSO force field expressions and units.

        Example with an LJ expression

        'Residue_1': {
        'expression': '4 * epsilon * ((sigma/r)**12 - (sigma/r)**6)',
        'expression_form': 'LJ',
        'expression_scalar': 1.0,
        'q_units': 'coulomb',
        'sigma_units': 'nm',
        'epsilon_units': 'kJ/mol'}
        }
    """
    eqn_gomc_std_forms_dict = {
        "LJ": "4*epsilon * ((sigma/r)**12 - (sigma/r)**6)",
        "Mie": "(n/(n-m)) * (n/m)**(m/(n-m)) * epsilon * ((sigma/r)**n - (sigma/r)**m)",
        "Exp6": "epsilon*alpha/(alpha-6) * (6/alpha*exp(alpha*(1-r/Rmin)) - (Rmin/r)**6)",
    }

    atomtypes_data_expression_data_dict = {}
    for res_i in atom_types_dict.keys():
        for atom_type_m in atom_types_dict[res_i]["atom_types"]:
            modified_atom_type_iter = f"{res_i}_{atom_type_m.name}"
            atomtypes_data_dict_iter = {
                modified_atom_type_iter: {
                    "expression": None,
                    "expression_form": None,
                    "expression_scalar": None,
                }
            }
            expression_iter = atom_types_dict[res_i]["expression"]
            if (
                atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_form"
                ]
                is None
                and atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_scalar"
                ]
                is None
            ):
                [
                    form_output,
                    form_scalar,
                ] = evaluate_nonbonded_lj_format_with_scaler(
                    expression_iter, eqn_gomc_std_forms_dict["LJ"]
                )
                atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression"
                ] = expression_iter
                atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_form"
                ] = form_output
                atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_scalar"
                ] = form_scalar

            if (
                atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_form"
                ]
                is None
                and atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_scalar"
                ]
                is None
            ):
                [
                    form_output,
                    form_scalar,
                ] = evaluate_nonbonded_mie_format_with_scaler(
                    expression_iter, eqn_gomc_std_forms_dict["Mie"]
                )
                atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression"
                ] = expression_iter
                atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_form"
                ] = form_output
                atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_scalar"
                ] = form_scalar

            if (
                atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_form"
                ]
                is None
                and atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_scalar"
                ]
                is None
            ):
                [
                    form_output,
                    form_scalar,
                ] = evaluate_nonbonded_exp6_format_with_scaler(
                    expression_iter, eqn_gomc_std_forms_dict["Exp6"]
                )
                atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression"
                ] = expression_iter
                atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_form"
                ] = form_output
                atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_scalar"
                ] = form_scalar

            if (
                atomtypes_data_dict_iter[modified_atom_type_iter]["expression"]
                is None
                or atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_form"
                ]
                is None
                or atomtypes_data_dict_iter[modified_atom_type_iter][
                    "expression_scalar"
                ]
                is None
            ):
                print_error_text = (
                    "ERROR: the {} residue does not match the listed standard or scaled "
                    "LJ, Mie, or Exp6 expressions in the {} function."
                    "".format(
                        modified_atom_type_iter,
                        "get_atom_type_expressions_and_scalars",
                    )
                )
                raise ValueError(print_error_text)

            atomtypes_data_expression_data_dict.update(atomtypes_data_dict_iter)

    return atomtypes_data_expression_data_dict


# compare harmonic bond equations or expressions
def evaluate_harmonic_bond_format_with_scaler(new_bond_form, base_bond_form):
    """Compare a new harmonic bond form to a base harmonic bond form (new bond form / base bond form).

    If the new form is the same as the base form, other than a scaling factor,
    it labels the form and provides its scalar from the base form.

    Parameters
    ----------
    new_bond_form : str
        The new bond form that will be compared or divided by the base form.
    base_bond_form : str
        The base bond form, which is the standard form.

    Returns
    -------
    list, [form_output, form_scalar]
        form_output : str
            'HarmonicBondPotential', if the new_bond_form variable is a harmonic bond.
            None, if the new_bond_form variable is not a harmonic bond.
        form_scalar : float
            float, if the new_bond_form variable is a harmonic bond.
            None, if the new_bond_form variable is not a harmonic bond.
    """
    try:
        eqn_ratio, k, r, r_eq = sympy.symbols("eqn_ratio k r r_eq")
        values = sympy.nonlinsolve(
            [
                eqn_ratio
                - sympy.sympify(new_bond_form) / sympy.sympify(base_bond_form),
            ],
            [eqn_ratio],
        )

        form_scalar = float(list(values)[0][0])
        form_output = "HarmonicBondPotential"
    except:
        form_scalar = None
        form_output = None

    return [form_output, form_scalar]


# compare harmonic angle equations or expressions
def evaluate_harmonic_angle_format_with_scaler(new_angle_form, base_angle_form):
    """Compare a new harmonic angle form to a base harmonic angle form (new angle form / base angle form).

    If the new form is the same as the base form, other than a scaling factor,
    it labels the form and provides its scalar from the base form.

    Parameters
    ----------
    new_angle_form : str
        The new angle form that will be compared or divided by the base form.
    base_angle_form : str
        The base angle form, which is the standard form.

    Returns
    -------
    list, [form_output, form_scalar]
        form_output : str
            'HarmonicAnglePotential', if the new_angle_form variable a harmonic angle.
            None, if the new_angle_form variable is not a harmonic angle.
        form_scalar : float
            float, if the new_angle_form variable is a harmonic angle.
            None, if the new_angle_form variable is not a harmonic.
    """
    try:
        eqn_ratio, k, theta, theta_eq = sympy.symbols(
            "eqn_ratio k theta theta_eq"
        )
        values = sympy.nonlinsolve(
            [
                eqn_ratio
                - sympy.sympify(new_angle_form)
                / sympy.sympify(base_angle_form),
            ],
            [eqn_ratio],
        )

        form_scalar = float(list(values)[0][0])
        form_output = "HarmonicAnglePotential"

    except:
        form_scalar = None
        form_output = None

    return [form_output, form_scalar]


# check for the harmonic torsion potential equations or expressions
def evaluate_harmonic_torsion_format_with_scaler(
    new_torsion_form, base_torsion_form
):
    """Compare a new harmonic torsion form to a base harmonic torsion form (new torsion form / base torsion form).

    If the new form is the same as the base form, other than a scaling factor,
    it labels the form and provides its scalar from the base form.

    Parameters
    ----------
    new_torsion_form : str
        The new harmonic torsion form that will be compared or divided by the base form.
    base_torsion_form : str
        The base harmonic torsion form, which is the standard form.

    Returns
    -------
    list, [form_output, form_scalar]
        form_output : str
            'HarmonicTorsionPotential', if the new_torsion_form variable is a harmonic torsion.
            None, if the new_torsion_form variable is not a harmonic torsion.
        form_scalar : float
            float, if the new_torsion_form variable is a harmonic torsion.
            None, if the new_torsion_form variable is not a harmonic torsion.
    """
    try:
        eqn_ratio, k, phi, phi_eq = sympy.symbols("eqn_ratio k phi phi_eq")
        values = sympy.nonlinsolve(
            [
                eqn_ratio
                - sympy.sympify(new_torsion_form)
                / sympy.sympify(base_torsion_form),
            ],
            [eqn_ratio],
        )

        form_scalar = float(list(values)[0][0])
        form_output = "HarmonicTorsionPotential"

    except:
        form_scalar = None
        form_output = None

    return [form_output, form_scalar]


# check for the OPLS torsion potential equations or expressions
def evaluate_OPLS_torsion_format_with_scaler(
    new_torsion_form, base_torsion_form
):
    """Compare a new OPLS torsion form to a base OPLS torsion form (new torsion form / base torsion form).

    If the new form is the same as the base form, other than a scaling factor,
    it labels the form and provides its scalar from the base form.

    Parameters
    ----------
    new_torsion_form : str
        The new OPLS torsion form that will be compared or divided by the base form.
    base_torsion_form : str
        The base OPLS torsion form, which is the standard form.

    Returns
    -------
    list, [form_output, form_scalar]
        form_output : str
            'OPLSTorsionPotential', if the new_torsion_form variable is an OPLS torsion.
            None, if the new_torsion_form variable is not an OPLS torsion.
        form_scalar : float
            float, if the new_torsion_form variable is an OPLS torsion.
            None, if the new_torsion_form variable is not an OPLS torsion.
    """
    try:
        eqn_ratio, k0, k1, k2, k3, k4, phi = sympy.symbols(
            "eqn_ratio k0 k1 k2 k3 k4 phi"
        )
        values = sympy.nonlinsolve(
            [
                eqn_ratio
                - sympy.sympify(new_torsion_form)
                / sympy.sympify(base_torsion_form),
            ],
            [eqn_ratio],
        )

        form_scalar = float(list(values)[0][0])
        form_output = "OPLSTorsionPotential"

    except:
        form_scalar = None
        form_output = None

    return [form_output, form_scalar]


# check for the periodic torsion potential equations or expressions
def evaluate_periodic_torsion_format_with_scaler(
    new_torsion_form, base_torsion_form
):
    """Compare a new periodic torsion form to a base periodic torsion form (new torsion form / base torsion form).

    If the new form is the same as the base form, other than a scaling factor,
    it labels the form and provides its scalar from the base form.

    Parameters
    ----------
    new_torsion_form : str
        The new periodic torsion form that will be compared or divided by the base formv
    base_torsion_form : str
        The base periodic torsion form, which is the standard form.

    Returns
    -------
    list, [form_output, form_scalar]
        form_output : str
            'PeriodicTorsionPotential', if the new_torsion_form variable is a periodic torsion.
            None, if the new_torsion_form variable is not a periodic torsion.
        form_scalar : float
            float, if the new_torsion_form variable is a periodic torsion.
            None, if the new_torsion_form variable is not a periodic torsion.
    """
    try:
        eqn_ratio, k, n, phi, phi_eq = sympy.symbols("eqn_ratio k n phi phi_eq")
        values = sympy.nonlinsolve(
            [
                eqn_ratio
                - sympy.sympify(new_torsion_form)
                / sympy.sympify(base_torsion_form),
            ],
            [eqn_ratio],
        )

        form_scalar = float(list(values)[0][0])
        form_output = "PeriodicTorsionPotential"

    except:
        form_scalar = None
        form_output = None

    return [form_output, form_scalar]


# check for the RyckaertBellemans (RB) torsion potential equations or expressions
def evaluate_RB_torsion_format_with_scaler(new_torsion_form, base_torsion_form):
    """Compare a new Ryckaert-Bellemans (RB) torsion form to a base torsion form (new torsion form / base torsion form).

    If the new form is the same as the base form, other than a scaling factor,
    it labels the form and provides its scalar from the base form.

    Parameters
    ----------
    new_torsion_form : str
        The new RB torsion form that will be compared or divided by the base formv
    base_torsion_form : str
        The base RB torsion form, which is the standard form.

    Returns
    -------
    list, [form_output, form_scalar]
        form_output : str
            'RyckaertBellemansTorsionPotential', if the new_torsion_form variable is an RB torsion.
            None, if the new_torsion_form variable is not an RB torsion.
        form_scalar : float
            float, if the new_torsion_form variable is an RB torsion.
            None, if the new_torsion_form variable is not an RB torsion.
    """
    try:
        eqn_ratio, c0, c1, c2, c3, c4, c5, psi = sympy.symbols(
            "eqn_ratio c0 c1 c2 c3 c4 c5 psi"
        )
        values = sympy.nonlinsolve(
            [
                eqn_ratio
                - sympy.sympify(new_torsion_form)
                / sympy.sympify(base_torsion_form),
            ],
            [eqn_ratio],
        )

        form_scalar = float(list(values)[0][0])
        form_output = "RyckaertBellemansTorsionPotential"

    except:
        form_scalar = None
        form_output = None

    return [form_output, form_scalar]


# check for the harmonic improper potential equations or expressions
def evaluate_harmonic_improper_format_with_scaler(
    new_improper_form, base_improper_form
):
    """Compare a new harmonic improper form to a base harmonic improper form (new improper form / base improper form).

    If the new form is the same as the base form, other than a scaling factor,
    it labels the form and provides its scalar from the base form.

    Parameters
    ----------
    new_improper_form : str
        The new harmonic improper form that will be compared or divided by the base form.
    base_improper_form : str
        The base harmonic improper form, which is the standard form.

    Returns
    -------
    list, [form_output, form_scalar]
        form_output : str
            'HarmonicImproperPotential', if the new_improper_form variable is a harmonic improper.
            None, if the new_improper_form variable is not a harmonic improper.
        form_scalar : float
            float, if the new_improper_form variable is a harmonic improper.
            None, if the new_improper_form variable is not a harmonic improper.
    """
    try:
        eqn_ratio, k, phi, phi_eq = sympy.symbols("eqn_ratio k phi phi_eq")
        values = sympy.nonlinsolve(
            [
                eqn_ratio
                - sympy.sympify(new_improper_form)
                / sympy.sympify(base_improper_form),
            ],
            [eqn_ratio],
        )

        form_scalar = float(list(values)[0][0])
        form_output = "HarmonicImproperPotential"

    except:
        form_scalar = None
        form_output = None

    return [form_output, form_scalar]


# check for the periodic improper potential equations or expressions
def evaluate_periodic_improper_format_with_scaler(
    new_improper_form, base_improper_form
):
    """Compare a new periodic improper form to a base periodic improper form (new improper form / base improper form).

    If the new form is the same as the base form, other than a scaling factor,
    it labels the form and provides its scalar from the base form.

    Parameters
    ----------
    new_improper_form : str
        The new periodic improper form that will be compared or divided by the base form.
    base_improper_form : str
        The base periodic improper form, which is the standard form.

    Returns
    -------
    list, [form_output, form_scalar]
        form_output : str
            'PeriodicImproperPotential', if the new_improper_form variable is a periodic improper.
            None, if the new_improper_form variable is not a periodic improper.
        form_scalar : float
            float, if the new_improper_form variable is a periodic improper.
            None, if the new_improper_form variable is not a periodic improper.
    """
    try:
        eqn_ratio, k, n, phi, phi_eq = sympy.symbols("eqn_ratio k n phi phi_eq")
        values = sympy.nonlinsolve(
            [
                eqn_ratio
                - sympy.sympify(new_improper_form)
                / sympy.sympify(base_improper_form),
            ],
            [eqn_ratio],
        )

        form_scalar = float(list(values)[0][0])
        form_output = "PeriodicImproperPotential"

    except:
        form_scalar = None
        form_output = None

    return [form_output, form_scalar]
