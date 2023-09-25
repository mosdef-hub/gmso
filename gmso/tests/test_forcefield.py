import lxml
import pytest
import unyt as u
from lxml.etree import DocumentInvalid
from sympy import sympify
from unyt.testing import assert_allclose_units

from gmso.core.forcefield import ForceField
from gmso.core.improper_type import ImproperType
from gmso.exceptions import (
    ForceFieldError,
    ForceFieldParseError,
    GMSOError,
    MissingAtomTypesError,
    MissingPotentialError,
)
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import allclose_units_mixed, get_path


class TestForceField(BaseTest):
    @pytest.fixture
    def ff(self):
        return ForceField(get_path("ff-example0.xml"))

    @pytest.fixture
    def named_groups_ff(self):
        return ForceField(get_path("ff-example1.xml"))

    @pytest.fixture
    def opls_ethane_foyer(self):
        return ForceField(
            get_path(filename=get_path("oplsaa-ethane_foyer.xml"))
        )

    @pytest.fixture(scope="session")
    def non_element_ff(self):
        return ForceField(get_path(filename="non-element-type-ff.xml"))

    def test_ff_name_version_from_xml(self, ff):
        assert ff.name == "ForceFieldOne"
        assert ff.version == "0.4.1"

    def test_scaling_factors_from_xml(self, ff):
        assert ff.scaling_factors["nonBonded14Scale"] == 0.67
        assert ff.scaling_factors["electrostatics14Scale"] == 0.5

    def test_ff_combining_rule(self, ff, opls_ethane_foyer):
        assert ff.combining_rule == "lorentz"
        assert opls_ethane_foyer.combining_rule == "geometric"

    @pytest.mark.parametrize(
        "unit_name,unit_value",
        [
            ("energy", u.Unit(u.kb)),
            ("mass", u.gram / u.mol),
            ("charge", u.coulomb),
            ("distance", u.nm),
        ],
    )
    def test_units_from_xml(self, ff, unit_name, unit_value):
        assert len(ff.units.keys()) == 4
        assert ff.units[unit_name] == unit_value

    def test_ff_atomtypes_from_xml(self, ff):
        assert len(ff.atom_types) == 3
        assert "Ar" in ff.atom_types
        assert "Xe" in ff.atom_types
        assert ff.atom_types["Ar"].get_tag("element") == "Ar"
        assert ff.atom_types["Xe"].get_tag("element") == "Xe"

        assert sympify("r") in ff.atom_types["Ar"].independent_variables
        assert ff.atom_types["Ar"].parameters["A"] == u.unyt_quantity(
            0.1, u.kcal / u.mol
        )
        assert ff.atom_types["Ar"].parameters["B"] == u.unyt_quantity(4.0, u.nm)
        assert ff.atom_types["Ar"].parameters["C"] == u.unyt_quantity(
            0.5, u.kcal / u.mol * u.nm**6
        )
        assert ff.atom_types["Ar"].mass == u.unyt_quantity(39.948, u.amu)
        assert ff.atom_types["Ar"].charge == u.unyt_quantity(0.0, u.coulomb)
        assert ff.atom_types["Ar"].description == "Argon atom"
        assert ff.atom_types["Ar"].definition == "Ar"
        assert ff.atom_types["Ar"].expression == sympify(
            "(A*exp(-B/r) - C/r**6)"
        )

        assert sympify("r") in ff.atom_types["Xe"].independent_variables
        assert "A" in ff.atom_types["Xe"].parameters
        assert ff.atom_types["Xe"].parameters["A"] == u.unyt_quantity(
            0.2, u.kcal / u.mol
        )
        assert ff.atom_types["Xe"].parameters["B"] == u.unyt_quantity(5.0, u.nm)
        assert ff.atom_types["Xe"].parameters["C"] == u.unyt_quantity(
            0.3, u.kcal / u.mol * u.nm**6
        )
        assert ff.atom_types["Xe"].mass == u.unyt_quantity(131.293, u.amu)
        assert ff.atom_types["Xe"].charge == u.unyt_quantity(0.0, u.coulomb)
        assert ff.atom_types["Xe"].description == "Xenon atom"
        assert ff.atom_types["Xe"].definition == "Xe"
        assert ff.atom_types["Xe"].expression == sympify(
            "(A*exp(-B/r) - C/r**6)"
        )

        assert ff.atom_types["Li"].charge == u.unyt_quantity(1.0, u.coulomb)

    def test_ff_bondtypes_from_xml(self, ff):
        assert len(ff.bond_types) == 2
        assert "Ar~Ar" in ff.bond_types
        assert "Xe~Xe" in ff.bond_types

        assert sympify("r") in ff.bond_types["Ar~Ar"].independent_variables
        assert ff.bond_types["Ar~Ar"].parameters["r_eq"] == u.unyt_quantity(
            10.0, u.nm
        )
        assert ff.bond_types["Ar~Ar"].parameters["k"] == u.unyt_quantity(
            10000, u.kJ / u.mol
        )
        assert ff.bond_types["Ar~Ar"].member_types == ("Ar", "Ar")

        assert sympify("r") in ff.bond_types["Xe~Xe"].independent_variables
        assert ff.bond_types["Xe~Xe"].parameters["r_eq"] == u.unyt_quantity(
            10.0, u.nm
        )
        assert ff.bond_types["Xe~Xe"].parameters["k"] == u.unyt_quantity(
            20000, u.kJ / u.mol
        )
        assert ff.bond_types["Xe~Xe"].member_types == ("Xe", "Xe")

    def test_ff_angletypes_from_xml(self, ff):
        assert len(ff.angle_types) == 2
        assert "Ar~Ar~Ar" in ff.angle_types
        assert "Xe~Xe~Xe" in ff.angle_types

        assert sympify("r") in ff.angle_types["Ar~Ar~Ar"].independent_variables
        assert ff.angle_types["Ar~Ar~Ar"].parameters["r_eq"] == u.unyt_quantity(
            10.0, u.nm
        )
        assert ff.angle_types["Ar~Ar~Ar"].parameters["z"] == u.unyt_quantity(
            100, u.kJ / u.mol
        )
        assert ff.angle_types["Ar~Ar~Ar"].member_types == ("Ar", "Ar", "Ar")

        assert sympify("r") in ff.angle_types["Xe~Xe~Xe"].independent_variables
        assert ff.angle_types["Xe~Xe~Xe"].parameters["r_eq"] == u.unyt_quantity(
            10.0, u.nm
        )
        assert ff.angle_types["Xe~Xe~Xe"].parameters["z"] == u.unyt_quantity(
            20, u.kJ / u.mol
        )
        assert ff.angle_types["Xe~Xe~Xe"].member_classes == ("Xe", "Xe", "Xe")

    def test_ff_dihedraltypes_from_xml(self, ff):
        assert len(ff.dihedral_types) == 2
        assert "Xe~Xe~Xe~Xe" in ff.dihedral_types
        assert "Ar~Ar~Ar~Ar" in ff.dihedral_types

        assert (
            sympify("r")
            in ff.dihedral_types["Ar~Ar~Ar~Ar"].independent_variables
        )
        assert ff.dihedral_types["Ar~Ar~Ar~Ar"].parameters[
            "r_eq"
        ] == u.unyt_quantity(10.0, u.nm)
        assert ff.dihedral_types["Ar~Ar~Ar~Ar"].parameters[
            "z"
        ] == u.unyt_quantity(100, u.kJ / u.mol)
        assert ff.dihedral_types["Ar~Ar~Ar~Ar"].member_classes == (
            "Ar",
            "Ar",
            "Ar",
            "Ar",
        )

        assert (
            sympify("r")
            in ff.dihedral_types["Xe~Xe~Xe~Xe"].independent_variables
        )
        assert ff.dihedral_types["Xe~Xe~Xe~Xe"].parameters[
            "r_eq"
        ] == u.unyt_quantity(10.0, u.nm)
        assert ff.dihedral_types["Xe~Xe~Xe~Xe"].parameters[
            "z"
        ] == u.unyt_quantity(20, u.kJ / u.mol)
        assert ff.dihedral_types["Xe~Xe~Xe~Xe"].member_classes == (
            "Xe",
            "Xe",
            "Xe",
            "Xe",
        )

    def test_ff_impropertypes_from_xml(self, ff):
        assert len(ff.improper_types) == 1
        assert "Xe~Xe~Xe~Xe" in ff.improper_types

        assert (
            sympify("r")
            in ff.improper_types["Xe~Xe~Xe~Xe"].independent_variables
        )
        assert ff.improper_types["Xe~Xe~Xe~Xe"].parameters[
            "r_eq"
        ] == u.unyt_quantity(10.0, u.nm)
        assert ff.improper_types["Xe~Xe~Xe~Xe"].parameters[
            "z"
        ] == u.unyt_quantity(20, u.kJ / u.mol)
        assert ff.improper_types["Xe~Xe~Xe~Xe"].member_types == (
            "Xe",
            "Xe",
            "Xe",
            "Xe",
        )

    def test_ff_pairpotentialtypes_from_xml(self, ff):
        assert len(ff.pairpotential_types) == 1
        assert "Xe~Xe" in ff.pairpotential_types

        assert (
            sympify("r")
            in ff.pairpotential_types["Xe~Xe"].independent_variables
        )
        assert ff.pairpotential_types["Xe~Xe"].parameters[
            "sigma"
        ] == u.unyt_quantity(10.0, u.nm)
        assert ff.pairpotential_types["Xe~Xe"].parameters[
            "k"
        ] == u.unyt_quantity(0.1, u.kJ / u.mol)
        assert ff.pairpotential_types["Xe~Xe"].member_types == ("Xe", "Xe")

    def test_ff_charmm_xml(self):
        charm_ff = ForceField(get_path("trimmed_charmm.xml"), backend="ffutils")

        assert charm_ff.name == "topologyCharmm"
        assert "*~CS~SS~*" in charm_ff.dihedral_types

        assert isinstance(
            charm_ff.dihedral_types["*~CE1~CE1~*"].parameters["k"], u.unyt_array
        )

        # This ensures that even though the parameters is a list, they can be hashed (by equality checks)
        assert (
            charm_ff.dihedral_types["*~CE1~CE1~*"]
            == charm_ff.dihedral_types["*~CE1~CE1~*"]
        )
        assert len(charm_ff.dihedral_types["*~CE1~CE1~*"].parameters["k"]) == 2

        # Test Correct Parameter Values
        assert_allclose_units(
            charm_ff.dihedral_types["*~CE1~CE1~*"].parameters["k"],
            [0.6276, 35.564] * u.kJ,
            rtol=1e-5,
            atol=1e-8,
        )

    def test_non_unique_params(self):
        # TODO: this should throw this error from forcefield-utilties, but currently does not.
        # with pytest.raises(DocumentInvalid):
        #    ForceField(get_path("ff-example-nonunique-params.xml"))
        pass

    def test_missing_params(self):
        # TODO: raise same error if backend loader is forcefield-utilities
        with pytest.raises(ForceFieldParseError):
            ForceField(
                get_path("ff-example-missing-parameter.xml"), backend="gmso"
            )

    def test_elementary_charge_to_coulomb(self, ff):
        elementary_charge = ff.atom_types["Li"].charge.to(u.elementary_charge)
        assert elementary_charge.units == u.Unit(u.elementary_charge)

    def test_atomclass_groups_charm_buck_ff(self):
        ff = ForceField(get_path("opls_charmm_buck.xml"))
        assert len(ff.atom_class_groups["CT"]) == 2

    def test_ff_periodic_dihedrals_from_alphanumeric_symbols(self):
        ff = ForceField(get_path("opls_charmm_buck.xml"))
        assert "A" in ff.atom_types["buck_O"].parameters
        with pytest.raises(TypeError):
            assert len(
                ff.dihedral_types["opls_140~*~*~opls_140"].parameters["c0"]
            )
        assert ff.dihedral_types["NH2~CT1~C~O"].parameters[
            "delta"
        ] == u.unyt_quantity(0.0, "degree")

    def test_ff_from_etree(self):
        # TODO: load using backend forcefield-utilities from etree
        ff_etree = lxml.etree.parse(get_path("opls_charmm_buck.xml"))
        ff = ForceField(ff_etree, backend="gmso")
        assert ff

    def test_ff_from_etree_iterable(self):
        # TODO: load using backend forcefield-utilities from etree
        ff_etrees = [
            lxml.etree.parse(get_path("opls_charmm_buck.xml")),
            lxml.etree.parse(get_path("trimmed_charmm.xml")),
        ]
        ff = ForceField(ff_etrees, backend="gmso")
        assert ff

    def test_ff_mixed_type_error(self):
        with pytest.raises(TypeError):
            ff = ForceField([5, "20"])

    def test_named_potential_groups(self):
        # TODO: get potential groups using backend forcefield-utilities
        named_groups_ff = ForceField(
            get_path("ff-example1.xml"), backend="gmso"
        )
        assert named_groups_ff.potential_groups["BuckinghamPotential"]
        assert (
            named_groups_ff.angle_types["Xe~Xe~Xe"]
            in named_groups_ff.potential_groups["HarmonicAngle"].values()
        )
        assert len(named_groups_ff.potential_groups["BuckinghamPotential"]) == 3
        assert len(named_groups_ff.potential_groups["HarmonicBond"]) == 2
        assert len(named_groups_ff.potential_groups["HarmonicAngle"]) == 2
        assert len(named_groups_ff.potential_groups["PeriodicProper"]) == 2
        assert len(named_groups_ff.potential_groups["RBProper"]) == 1
        assert len(named_groups_ff.potential_groups["LJ"]) == 1

    def test_potential_types_by_expression(self, named_groups_ff):
        atom_types_grouped_by_expression = (
            named_groups_ff.group_atom_types_by_expression()
        )
        bond_types_grouped_by_expression = (
            named_groups_ff.group_bond_types_by_expression()
        )
        angle_types_grouped_by_expression = (
            named_groups_ff.group_angle_types_by_expression()
        )
        dihedral_types_grouped_by_expression = (
            named_groups_ff.group_dihedral_types_by_expression()
        )
        improper_types_grouped_by_expression = (
            named_groups_ff.group_improper_types_by_expression()
        )
        pairpotential_types_grouped_by_expression = (
            named_groups_ff.group_pairpotential_types_by_expression()
        )

        assert (
            len(atom_types_grouped_by_expression["A*exp(-B/r) - C/r**6"]) == 3
        )
        assert len(bond_types_grouped_by_expression["0.5*k*(r - r_eq)**2"]) == 2
        assert (
            len(angle_types_grouped_by_expression["0.5*z*(r - r_eq)**2"]) == 2
        )
        assert (
            len(dihedral_types_grouped_by_expression["0.5*z*(r - r_eq)**2"])
            == 2
        )
        assert (
            len(improper_types_grouped_by_expression["0.5*z*(r - r_eq)**2"])
            == 1
        )
        assert (
            len(
                pairpotential_types_grouped_by_expression[
                    "4*k*(-sigma**6/r**6 + sigma**12/r**12)"
                ]
            )
            == 1
        )

    def test_forcefield_missing_atom_types(self):
        with pytest.raises(MissingAtomTypesError):
            ff = ForceField(
                get_path(filename=get_path("ff_missing_atom_types.xml")),
                backend="gmso",
            )

    def test_forcefield_missing_atom_types_non_strict(self):
        ff = ForceField(
            get_path(filename=get_path("ff_missing_atom_types.xml")),
            strict=False,
            backend="gmso",
        )

    def test_forcefeld_get_potential_atom_type(self, opls_ethane_foyer):
        at = opls_ethane_foyer.get_potential("atom_type", key=["opls_135"])
        assert at.expression == sympify(
            "ep * ((sigma/r)**12 - (sigma/r)**6) + q / (e0 * r)"
        )

        params = at.parameters
        assert "ep" in params
        assert "sigma" in params
        assert "e0" in params
        assert sympify("r") in at.independent_variables

        assert allclose_units_mixed(
            params.values(),
            [
                0.276144 * u.kJ / u.mol,
                0.35 * u.nm,
                8.8542e-12 * u.Unit("A**2*s**4/(kg*m**3)"),
                -0.18 * u.C,
            ],
        )

    def test_forcefield_get_parameters_atom_type(self, opls_ethane_foyer):
        params = opls_ethane_foyer.get_parameters("atom_type", key=["opls_140"])

        assert allclose_units_mixed(
            params.values(),
            [
                0.12552 * u.kJ / u.mol,
                0.25 * u.nm,
                8.8542e-12 * u.Unit("A**2*s**4/(kg*m**3)"),
                0.06 * u.C,
            ],
        )

    def test_forcefield_get_parameters_atom_type_copy(self, opls_ethane_foyer):
        params = opls_ethane_foyer.get_parameters(
            "atom_type", key=["opls_140"], copy=False
        )
        params_copy = opls_ethane_foyer.get_parameters(
            "atom_type", key=["opls_140"], copy=True
        )
        assert allclose_units_mixed(params.values(), params_copy.values())

    def test_forcefield_get_potential_bond_type(self, opls_ethane_foyer):
        bt = opls_ethane_foyer.get_potential(
            "bond_type", key=["opls_135", "opls_140"]
        )
        assert bt.name == "BondType-Harmonic-2"
        params = bt.parameters
        assert "k" in params
        assert "r_eq" in params

        assert sympify("r") in bt.independent_variables

        assert allclose_units_mixed(
            params.values(), [284512.0 * u.kJ / u.nm**2, 0.109 * u.nm]
        )

    def test_forcefield_get_potential_bond_type_reversed(
        self, opls_ethane_foyer
    ):
        assert opls_ethane_foyer.get_potential(
            "bond_type", ["opls_135", "opls_140"]
        ) == opls_ethane_foyer.get_potential(
            "bond_type", ["opls_140", "opls_135"]
        )

    def test_forcefield_get_parameters_bond_type(self, opls_ethane_foyer):
        params = opls_ethane_foyer.get_parameters(
            "bond_type", key=["opls_135", "opls_135"]
        )

        assert allclose_units_mixed(
            params.values(), [224262.4 * u.kJ / u.nm**2, 0.1529 * u.nm]
        )

    def test_forcefield_get_potential_angle_type(self, opls_ethane_foyer):
        at = opls_ethane_foyer.get_potential(
            "angle_type", key=["opls_135", "opls_135", "opls_140"]
        )
        assert at.name == "AngleType-Harmonic-1"
        params = at.parameters
        assert "k" in params
        assert "theta_eq" in params

        assert sympify("theta") in at.independent_variables

        assert allclose_units_mixed(
            params.values(),
            [313.8 * u.kJ / u.radian**2, 1.932079482 * u.radian],
        )

    def test_forcefield_get_potential_angle_type_reversed(
        self, opls_ethane_foyer
    ):
        assert opls_ethane_foyer.get_potential(
            "angle_type", ["opls_135", "opls_135", "opls_140"]
        ) == opls_ethane_foyer.get_potential(
            "angle_type", ["opls_140", "opls_135", "opls_135"]
        )

    def test_forcefield_get_parameters_angle_type(self, opls_ethane_foyer):
        params = opls_ethane_foyer.get_parameters(
            "angle_type", key=["opls_140", "opls_135", "opls_140"]
        )

        assert allclose_units_mixed(
            params.values(),
            [276.144 * u.kJ / u.radian**2, 1.8814649337 * u.radian],
        )

    def test_forcefield_get_potential_dihedral_type(self, opls_ethane_foyer):
        dt = opls_ethane_foyer.get_potential(
            "dihedral_type",
            key=["opls_140", "opls_135", "opls_135", "opls_140"],
        )
        assert dt.name == "DihedralType-RB-Proper-1"
        params = dt.parameters
        assert "c0" in params
        assert "c1" in params
        assert "c2" in params
        assert "c3" in params
        assert "c4" in params
        assert "c5" in params

        assert sympify("phi") in dt.independent_variables

        assert allclose_units_mixed(
            params.values(),
            [0.6276, 1.8828, 0.0, -2.5104, 0.0, 0.0] * u.kJ / u.mol,
        )

    def test_forcefield_get_parameters_dihedral_type(self, opls_ethane_foyer):
        params = opls_ethane_foyer.get_parameters(
            "dihedral_type",
            key=["opls_140", "opls_135", "opls_135", "opls_140"],
        )

        assert allclose_units_mixed(
            params.values(),
            [0.6276, 1.8828, 0.0, -2.5104, 0.0, 0.0] * u.kJ / u.mol,
        )

    def test_forcefield_get_potential_non_exisistent_group(
        self, opls_ethane_foyer
    ):
        with pytest.raises(ValueError):
            opls_ethane_foyer.get_potential("non_group", ["a", "b", "c"])

    def test_forcefield_get_potential_non_string_key(self, opls_ethane_foyer):
        with pytest.raises(TypeError):
            opls_ethane_foyer.get_potential("atom_type", key=[111])

    def test_get_atom_type_missing(self, opls_ethane_foyer):
        with pytest.raises(MissingPotentialError):
            opls_ethane_foyer._get_atom_type("opls_359", warn=False)

        with pytest.warns(UserWarning):
            opls_ethane_foyer._get_atom_type("opls_359", warn=True)

    def test_get_bond_type_missing(self, opls_ethane_foyer):
        with pytest.raises(MissingPotentialError):
            opls_ethane_foyer._get_bond_type(
                ["opls_359", "opls_600"], warn=False
            )

        with pytest.warns(UserWarning):
            opls_ethane_foyer._get_bond_type(
                ["opls_359", "opls_600"], warn=True
            )

    def test_get_angle_type_missing(self, opls_ethane_foyer):
        with pytest.raises(MissingPotentialError):
            opls_ethane_foyer._get_angle_type(
                ["opls_359", "opls_600", "opls_700"], warn=False
            )

        with pytest.warns(UserWarning):
            opls_ethane_foyer._get_angle_type(
                ["opls_359", "opls_600", "opls_700"], warn=True
            )

    def test_get_dihedral_type_missing(self, opls_ethane_foyer):
        with pytest.raises(MissingPotentialError):
            opls_ethane_foyer._get_dihedral_type(
                ["opls_359", "opls_600", "opls_700", "opls_800"], warn=False
            )

        with pytest.warns(UserWarning):
            opls_ethane_foyer._get_dihedral_type(
                ["opls_359", "opls_600", "opls_700", "opls_800"], warn=True
            )

    def test_get_improper_type_missing(self, opls_ethane_foyer):
        with pytest.raises(MissingPotentialError):
            opls_ethane_foyer._get_improper_type(
                ["opls_359", "opls_600", "opls_700", "opls_800"], warn=False
            )

        with pytest.warns(UserWarning):
            opls_ethane_foyer._get_improper_type(
                ["opls_359", "opls_600", "opls_700", "opls_800"], warn=True
            )

    def test_non_element_types(self, non_element_ff, opls_ethane_foyer):
        assert "_CH3" in non_element_ff.non_element_types
        assert "_CH2" in non_element_ff.non_element_types
        assert opls_ethane_foyer.non_element_types == set()
        assert len(opls_ethane_foyer.atom_types) > 0

        assert (
            non_element_ff.get_potential(
                group="atom_type", key="CH2_sp3"
            ).charge
            == 0
        )
        assert (
            non_element_ff.get_potential(
                group="atom_type", key="CH3_sp3"
            ).charge
            == 0
        )

        assert (
            non_element_ff.get_potential(
                group="atom_type", key="CH3_sp3"
            ).definition
            == "[_CH3;X1][_CH3,_CH2]"
        )
        assert (
            non_element_ff.get_potential(
                group="atom_type", key="CH2_sp3"
            ).definition
            == "[_CH2;X2]([_CH3,_CH2])[_CH3,_CH2]"
        )

    def test_forcefield_get_impropers_combinations(self):
        ff_with_impropers = ForceField()
        ff_with_impropers.name = "imp_ff"
        ff_with_impropers.improper_types = {
            "CT~CT~HC~HC": ImproperType(name="imp1"),
            "CT~HC~HC~HC": ImproperType(name="imp2"),
        }
        imp1 = ff_with_impropers.get_potential(
            "improper_type", ["CT", "HC", "HC", "CT"]
        )
        imp2 = ff_with_impropers.get_potential(
            "improper_type", ["CT", "HC", "CT", "HC"]
        )
        assert imp1.name == imp2.name
        assert imp1 is imp2

    def test_write_xml(self, opls_ethane_foyer):
        opls_ethane_foyer.to_xml("test_xml_writer.xml")
        reloaded_xml = ForceField("test_xml_writer.xml")
        get_names = lambda ff, param: [
            typed for typed in getattr(ff, param).keys()
        ]
        for param in [
            "atom_types",
            "bond_types",
            "angle_types",
            "dihedral_types",
        ]:
            assert get_names(opls_ethane_foyer, param) == get_names(
                reloaded_xml, param
            )

    def test_write_not_xml(self, opls_ethane_foyer):
        with pytest.raises(ForceFieldError):
            opls_ethane_foyer.to_xml("bad_path")

    def test_valid_sequence(self):
        for j in range(10):
            ff = ForceField(get_path("sequence_of_parameters_ff.xml"), "r")
            dih_with_list = ff.dihedral_types["*~C~C~*"]
            params = dih_with_list.get_parameters()
            assert u.allclose_units(params["theta_0"], [25, 32] * u.radian)
            assert u.allclose_units(params["k"], [38, 45] * u.kJ / u.mol)

    def test_deprecated_gmso(self):
        with pytest.warns(DeprecationWarning):
            ForceField(get_path("ff-example0.xml"), backend="gmso")

    def test_not_supoprted_backend(self, opls_ethane_foyer):
        # Unsupported ff parser backend
        with pytest.raises(GMSOError):
            ForceField(get_path("ff-example0.xml"), backend="bogus")

        # Unsupported ff writer backend
        with pytest.raises(NotImplementedError):
            opls_ethane_foyer.to_xml("test_xml_writer.xml", backend="ffutils")

        with pytest.raises(GMSOError):
            opls_ethane_foyer.to_xml("test_xml_writer.xml", backend="bogus")
