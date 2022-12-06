import mbuild as mb
import pytest
from mbuild.utils.io import has_foyer

from gmso.tests.base_test import BaseTest
from gmso.utils.equation_compare import (
    evaluate_harmonic_angle_format_with_scaler,
    evaluate_harmonic_bond_format_with_scaler,
    evaluate_harmonic_improper_format_with_scaler,
    evaluate_harmonic_torsion_format_with_scaler,
    evaluate_nonbonded_exp6_format_with_scaler,
    evaluate_nonbonded_lj_format_with_scaler,
    evaluate_nonbonded_mie_format_with_scaler,
    evaluate_OPLS_torsion_format_with_scaler,
    evaluate_periodic_improper_format_with_scaler,
    evaluate_periodic_torsion_format_with_scaler,
    evaluate_RB_torsion_format_with_scaler,
    get_atom_type_expressions_and_scalars,
)
from gmso.utils.io import get_fn
from gmso.utils.specific_ff_to_residue import specific_ff_to_residue

# base forms
input_base_lj_form = "4*epsilon * ((sigma/r)**12 - (sigma/r)**6)"
input_base_mie_form = (
    "(n/(n-m)) * (n/m)**(m/(n-m)) * epsilon * ((sigma/r)**n - (sigma/r)**m)"
)
input_base_exp6_form = (
    "epsilon*alpha/(alpha-6) * (6/alpha*exp(alpha*(1-r/Rmin)) - (Rmin/r)**6)"
)

# main_output_types
lj_output = "LJ"
mie_output = "Mie"
exp6_output = "Exp6"


@pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
class TestEqnCompare(BaseTest):
    def test_wrong_lj_equation_form(self):
        input_new_lj_form = "x"

        [form_output, form_scalar] = evaluate_nonbonded_lj_format_with_scaler(
            input_new_lj_form, input_base_lj_form
        )

        assert form_output == None
        assert form_scalar == None

    def test_scaled_lj_equation_form(self):
        input_new_lj_form = "8*epsilon * ((sigma/r)**12 - (sigma/r)**6)"

        [form_output, form_scalar] = evaluate_nonbonded_lj_format_with_scaler(
            input_new_lj_form, input_base_lj_form
        )

        assert form_output == lj_output
        assert form_scalar == 2.0

    def test_wrong_mie_equation_form(self):
        input_new_mie_form = "x"

        [form_output, form_scalar] = evaluate_nonbonded_mie_format_with_scaler(
            input_new_mie_form, input_base_mie_form
        )

        assert form_output == None
        assert form_scalar == None

    def test_scaled_mie_equation_form(self):
        input_new_mie_form = "2*(n/(n-m)) * (n/m)**(m/(n-m)) * epsilon * ((sigma/r)**n - (sigma/r)**m)"

        [form_output, form_scalar] = evaluate_nonbonded_mie_format_with_scaler(
            input_new_mie_form, input_base_mie_form
        )

        assert form_output == mie_output
        assert form_scalar == 2.0

    def test_wrong_exp6_equation_form(self):
        input_new_exp6_form = "x"

        [form_output, form_scalar] = evaluate_nonbonded_exp6_format_with_scaler(
            input_new_exp6_form, input_base_exp6_form
        )

        assert form_output == None
        assert form_scalar == None

    def test_scaled_exp6_equation_form(self):
        input_new_exp6_form = "2*epsilon*alpha/(alpha-6) * (6/alpha*exp(alpha*(1-r/Rmin)) - (Rmin/r)**6)"

        [form_output, form_scalar] = evaluate_nonbonded_exp6_format_with_scaler(
            input_new_exp6_form, input_base_exp6_form
        )

        assert form_output == exp6_output
        assert form_scalar == 2.0

    # the exp6 needs added when the standard exp6 is added to gmso overall.
    def test_find_lj_mie_exp6_forms_and_scalars(self):
        ethane_lj = mb.load(get_fn("ethane_ua.mol2"))
        ethane_lj.name = "ETH"
        ethane_mie = mb.load(get_fn("ethane_ua_mie.mol2"))
        ethane_mie.name = "ETHM"
        test_box = mb.fill_box(
            compound=[ethane_lj, ethane_mie], n_compounds=[1, 1], box=[4, 4, 4]
        )

        [
            test_topology,
            test_residues_applied_list,
            test_electrostatics14Scale_dict,
            test_nonBonded14Scale_dict,
            test_atom_types_dict,
            test_bond_types_dict,
            test_angle_types_dict,
            test_dihedral_types_dict,
            test_improper_types_dict,
            test_combining_rule_dict,
        ] = specific_ff_to_residue(
            test_box,
            forcefield_selection={
                "ETH": f"{get_fn('gmso_xmls/test_ffstyles/ethane_propane_ua_lorentz_combining.xml')}",
                "ETHM": f"{get_fn('gmso_xmls/test_ffstyles/ethane_propane_ua_Mie_lorentz_combining.xml')}",
            },
            residues=["ETH", "ETHM"],
            boxes_for_simulation=1,
        )

        atom_types_data_expression_data_dict = (
            get_atom_type_expressions_and_scalars(test_atom_types_dict)
        )

        assert (
            str(atom_types_data_expression_data_dict["ETH_CH3"]["expression"])
            == "4*epsilon*(-sigma**6/r**6 + sigma**12/r**12)"
        )
        assert (
            atom_types_data_expression_data_dict["ETH_CH3"]["expression_form"]
            == "LJ"
        )
        assert (
            atom_types_data_expression_data_dict["ETH_CH3"]["expression_scalar"]
            == 1.0
        )

        assert (
            str(atom_types_data_expression_data_dict["ETHM_CH3"]["expression"])
            == "epsilon*n*(n/m)**(m/(-m + n))*(-(sigma/r)**m + (sigma/r)**n)/(-m + n)"
        )
        assert (
            atom_types_data_expression_data_dict["ETHM_CH3"]["expression_form"]
            == "Mie"
        )
        assert (
            atom_types_data_expression_data_dict["ETHM_CH3"][
                "expression_scalar"
            ]
            == 1.0
        )

    # the exp6 needs added when the standard exp6 is added to gmso overall.
    def test_bad_eqn_for_find_lj_mie_exp6_forms_and_scalars(self):
        with pytest.raises(
            ValueError,
            match=r"ERROR: the ETHM_CH3 residue does not match the listed standard or scaled "
            r"LJ, Mie, or Exp6 expressions in the get_atom_type_expressions_and_scalars function.",
        ):
            ethane_lj = mb.load(get_fn("ethane_ua.mol2"))
            ethane_lj.name = "ETH"
            ethane_mie = mb.load(get_fn("ethane_ua_mie.mol2"))
            ethane_mie.name = "ETHM"
            test_box = mb.fill_box(
                compound=[ethane_lj, ethane_mie],
                n_compounds=[1, 1],
                box=[4, 4, 4],
            )

            [
                test_topology,
                test_residues_applied_list,
                test_electrostatics14Scale_dict,
                test_nonBonded14Scale_dict,
                test_atom_types_dict,
                test_bond_types_dict,
                test_angle_types_dict,
                test_dihedral_types_dict,
                test_improper_types_dict,
                test_combining_rule_dict,
            ] = specific_ff_to_residue(
                test_box,
                forcefield_selection={
                    "ETH": f"{get_fn('gmso_xmls/test_ffstyles/ethane_propane_ua_lorentz_combining.xml')}",
                    "ETHM": f"{get_fn('gmso_xmls/test_ffstyles/ethane_propane_ua_bad_eqn_lorentz_combining.xml')}",
                },
                residues=["ETH", "ETHM"],
                boxes_for_simulation=1,
            )

            atom_types_data_expression_data_dict = (
                get_atom_type_expressions_and_scalars(test_atom_types_dict)
            )

    # harmonic bond
    def test_wrong_harmonic_bond(self):
        input_base_harmonic_bond = "1 * k * (r-r_eq)**2"
        input_new_harmonic_form = "x"

        [form_output, form_scalar] = evaluate_harmonic_bond_format_with_scaler(
            input_new_harmonic_form, input_base_harmonic_bond
        )

        assert form_output == None
        assert form_scalar == None

    def test_scaled_harmonic_bond(self):
        # bond types
        input_base_harmonic_bond = "1 * k * (r-r_eq)**2"
        input_new_harmonic_form = "1/2 * k * (r-r_eq)**2"

        [form_output, form_scalar] = evaluate_harmonic_bond_format_with_scaler(
            input_new_harmonic_form, input_base_harmonic_bond
        )

        assert form_output == "HarmonicBondPotential"
        assert form_scalar == 0.5

    # harmonic angle
    def test_wrong_harmonic_angle(self):
        input_base_harmonic_angle = "1 * k * (theta - theta_eq)**2"
        input_new_harmonic_form = "x"

        [form_output, form_scalar] = evaluate_harmonic_angle_format_with_scaler(
            input_new_harmonic_form, input_base_harmonic_angle
        )

        assert form_output == None
        assert form_scalar == None

    def test_scaled_harmonic_angle(self):
        # angle types
        input_base_harmonic_angle = "1 * k * (theta - theta_eq)**2"
        input_new_harmonic_form = "1/2 * k * (theta - theta_eq)**2"

        [form_output, form_scalar] = evaluate_harmonic_angle_format_with_scaler(
            input_new_harmonic_form, input_base_harmonic_angle
        )

        assert form_output == "HarmonicAnglePotential"
        assert form_scalar == 0.5

    # harmonic torsion
    def test_wrong_harmonic_torsion(self):
        input_base_harmonic_torsion = "1 * k * (phi - phi_eq)**2"
        input_new_harmonic_form = "x"

        [
            form_output,
            form_scalar,
        ] = evaluate_harmonic_torsion_format_with_scaler(
            input_new_harmonic_form, input_base_harmonic_torsion
        )

        assert form_output == None
        assert form_scalar == None

    def test_scaled_harmonic_torsion(self):
        # torsion types
        input_base_harmonic_torsion = "1 * k * (phi - phi_eq)**2"
        input_new_harmonic_form = "1/2 * k * (phi - phi_eq)**2"

        [
            form_output,
            form_scalar,
        ] = evaluate_harmonic_torsion_format_with_scaler(
            input_new_harmonic_form, input_base_harmonic_torsion
        )

        assert form_output == "HarmonicTorsionPotential"
        assert form_scalar == 0.5

    # OPLS torsion
    def test_wrong_opls_torsion(self):
        input_base_opls_torsion = (
            "0.5 * k0 + "
            "0.5 * k1 * (1 + cos(phi)) + "
            "0.5 * k2 * (1 - cos(2*phi)) + "
            "0.5 * k3 * (1 + cos(3*phi)) + "
            "0.5 * k4 * (1 - cos(4*phi))"
        )
        input_new_opls_form = "x"

        [form_output, form_scalar] = evaluate_OPLS_torsion_format_with_scaler(
            input_new_opls_form, input_base_opls_torsion
        )

        assert form_output == None
        assert form_scalar == None

    def test_scaled_opls_torsion(self):
        # torsion types
        input_base_opls_torsion = (
            "1/2 * k0 + "
            "1/2 * k1 * (1 + cos(phi)) + "
            "1/2 * k2 * (1 - cos(2*phi)) + "
            "1/2 * k3 * (1 + cos(3*phi)) + "
            "1/2 * k4 * (1 - cos(4*phi))"
        )

        input_new_opls_form = (
            "2 * k0 + "
            "2 * k1 * (1 + cos(phi)) + "
            "2 * k2 * (1 - cos(2*phi)) + "
            "2 * k3 * (1 + cos(3*phi)) + "
            "2 * k4 * (1 - cos(4*phi))"
        )

        [form_output, form_scalar] = evaluate_OPLS_torsion_format_with_scaler(
            input_new_opls_form, input_base_opls_torsion
        )

        assert form_output == "OPLSTorsionPotential"
        assert form_scalar == 4

    # periodic torsion
    def test_wrong_periodic_torsion(self):
        input_base_periodic_torsion = "k * (1 + cos(n * phi - phi_eq))"
        input_new_periodic_form = "x"

        [
            form_output,
            form_scalar,
        ] = evaluate_periodic_torsion_format_with_scaler(
            input_new_periodic_form, input_base_periodic_torsion
        )

        assert form_output == None
        assert form_scalar == None

    def test_scaled_periodic_torsion(self):
        # torsion types
        input_base_periodic_torsion = "k * (1 + cos(n * phi - phi_eq))"

        input_new_periodic_form = "2 * k * (1 + cos(n * phi - phi_eq))"

        [
            form_output,
            form_scalar,
        ] = evaluate_periodic_torsion_format_with_scaler(
            input_new_periodic_form, input_base_periodic_torsion
        )

        assert form_output == "PeriodicTorsionPotential"
        assert form_scalar == 2

    # RB torsion
    def test_wrong_RB_torsion(self):
        input_base_RB_torsion = (
            "c0 * cos(phi)**0 + "
            "c1 * cos(phi)**1 + "
            "c2 * cos(phi)**2 + "
            "c3 * cos(phi)**3 + "
            "c4 * cos(phi)**4 + "
            "c5 * cos(phi)**5"
        )
        input_new_RB_form = "x"

        [form_output, form_scalar] = evaluate_RB_torsion_format_with_scaler(
            input_new_RB_form, input_base_RB_torsion
        )

        assert form_output == None
        assert form_scalar == None

    def test_scaled_RB_torsion(self):
        # torsion types
        input_base_RB_torsion = (
            "c0 * cos(phi)**0 + "
            "c1 * cos(phi)**1 + "
            "c2 * cos(phi)**2 + "
            "c3 * cos(phi)**3 + "
            "c4 * cos(phi)**4 + "
            "c5 * cos(phi)**5"
        )

        input_new_RB_form = (
            "2 * c0 * cos(phi)**0 + "
            "2 * c1 * cos(phi)**1 + "
            "2 * c2 * cos(phi)**2 + "
            "2 * c3 * cos(phi)**3 + "
            "2 * c4 * cos(phi)**4 + "
            "2 * c5 * cos(phi)**5"
        )

        [form_output, form_scalar] = evaluate_RB_torsion_format_with_scaler(
            input_new_RB_form, input_base_RB_torsion
        )

        assert form_output == "RyckaertBellemansTorsionPotential"
        assert form_scalar == 2

    # RB torsion
    def test_wrong_RB_torsion(self):
        input_base_harmonic_improper = "k * (phi - phi_eq)**2"

        input_new_harmonic_improper_form = "x"

        [
            form_output,
            form_scalar,
        ] = evaluate_harmonic_improper_format_with_scaler(
            input_new_harmonic_improper_form,
            input_base_harmonic_improper,
        )

        assert form_output == None
        assert form_scalar == None

    # scaled_harmonic_improper
    def test_harmonic_improper(self):
        input_base_harmonic_improper = "k * (phi - phi_eq)**2"
        input_new_harmonic_improper_form = "x"

        [
            form_output,
            form_scalar,
        ] = evaluate_harmonic_improper_format_with_scaler(
            input_new_harmonic_improper_form,
            input_base_harmonic_improper,
        )

        assert form_output == None
        assert form_scalar == None

    def test_scaled_harmonic_improper(self):
        # torsion types
        input_base_harmonic_improper = "k * (phi - phi_eq)**2"
        input_new_harmonic_improper_form = "2 * k * (phi - phi_eq)**2"

        [
            form_output,
            form_scalar,
        ] = evaluate_harmonic_improper_format_with_scaler(
            input_new_harmonic_improper_form,
            input_base_harmonic_improper,
        )

        assert form_output == "HarmonicImproperPotential"
        assert form_scalar == 2

    # scaled_periodic_improper
    def test_periodic_improper(self):
        input_base_periodic_improper = "k * (1 + cos(n * phi - phi_eq))"
        input_new_periodic_improper_form = "x"

        [
            form_output,
            form_scalar,
        ] = evaluate_periodic_improper_format_with_scaler(
            input_new_periodic_improper_form,
            input_base_periodic_improper,
        )

        assert form_output == None
        assert form_scalar == None

    def test_scaled_periodic_improper(self):
        # periodic types
        input_base_periodic_improper = "k * (1 + cos(n * phi - phi_eq))"
        input_new_periodic_improper_form = "2* k * (1 + cos(n * phi - phi_eq))"

        [
            form_output,
            form_scalar,
        ] = evaluate_periodic_improper_format_with_scaler(
            input_new_periodic_improper_form,
            input_base_periodic_improper,
        )

        assert form_output == "PeriodicImproperPotential"
        assert form_scalar == 2
