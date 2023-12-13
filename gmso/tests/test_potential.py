import pytest
import sympy
import unyt as u
from unyt.testing import assert_allclose_units

from gmso.core.atom import Atom
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType
from gmso.core.parametric_potential import ParametricPotential
from gmso.core.topology import Topology
from gmso.exceptions import GMSOError
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.tests.base_test import BaseTest
from gmso.utils.sorting import sort_by_classes, sort_by_types


class TestPotential(BaseTest):
    def test_new_potential(self):
        new_potential = ParametricPotential(
            name="mypotential",
            expression="a*x+b",
            parameters={"a": 1.0 * u.g, "b": 1.0 * u.m},
            independent_variables={"x"},
        )
        assert new_potential.name == "mypotential"
        assert new_potential.expression == sympy.sympify("a*x+b")
        assert_allclose_units(
            new_potential.parameters["a"], 1.0 * u.g, rtol=1e-5, atol=1e-8
        )
        assert_allclose_units(
            new_potential.parameters["b"], 1.0 * u.m, rtol=1e-5, atol=1e-8
        )
        assert new_potential.independent_variables == {sympy.symbols("x")}

    def test_setters(self):
        new_potential = ParametricPotential()
        new_potential.name = "SettingName"
        new_potential.set_expression(
            independent_variables=sympy.symbols({"r"}),
            expression="r*sigma*epsilon",
            parameters={
                "sigma": 1 * u.nm,
                "epsilon": 10 * u.Unit("kcal / mol"),
            },
        )

        assert new_potential.name == "SettingName"
        assert new_potential.independent_variables == {sympy.symbols("r")}
        assert new_potential.parameters == {
            "sigma": 1 * u.nm,
            "epsilon": 10 * u.Unit("kcal / mol"),
        }
        assert new_potential.expression == sympy.sympify("r * sigma * epsilon")

    def test_incorrect_indep_vars(self):
        with pytest.raises(ValueError):
            ParametricPotential(expression="x*y", independent_variables="z")

    def test_incorrect_expression(self):
        with pytest.raises(ValueError):
            ParametricPotential(name="mytype", expression=4.2)

    def test_expression_consistency(self):
        new_potential = ParametricPotential(
            name="mypotential",
            parameters={"x": 1 * u.m, "y": 10 * u.m},
            expression="x + y * z",
            independent_variables="z",
        )

        symbol_x, symbol_y, symbol_z = sympy.symbols("x y z")
        correct_expr = sympy.sympify("x+y*z")
        assert new_potential.expression.free_symbols == set(
            [symbol_x, symbol_y, symbol_z]
        )
        assert correct_expr == new_potential.expression

    def test_equivalance(self):
        ref_potential = ParametricPotential(
            name="mypotential",
            expression="a*x+b",
            parameters={"a": 1.0 * u.g, "b": 1.0 * u.m},
            independent_variables={"x"},
        )
        same_potential = ParametricPotential(
            name="mypotential",
            expression="a*x+b",
            parameters={"a": 1.0 * u.g, "b": 1.0 * u.m},
            independent_variables={"x"},
        )
        diff_name = ParametricPotential(
            name="nopotential",
            expression="a*x+b",
            parameters={"a": 1.0 * u.g, "b": 1.0 * u.m},
            independent_variables={"x"},
        )
        diff_expression = ParametricPotential(
            name="mypotential",
            expression="a+x*b",
            parameters={"a": 1.0 * u.g, "b": 1.0 * u.m},
            independent_variables={"x"},
        )
        diff_params = ParametricPotential(
            name="mypotential",
            expression="a*x+b",
            parameters={"a": 2.0 * u.g, "b": 1.0 * u.m},
            independent_variables={"x"},
        )

        assert ref_potential == same_potential
        assert ref_potential != diff_name
        assert ref_potential != diff_expression
        assert ref_potential != diff_params
        assert same_potential != diff_name
        assert same_potential != diff_expression
        assert same_potential != diff_params
        assert diff_name != diff_expression
        assert diff_expression != diff_params

    def test_set_expression(self):
        # Try changing the expression but keep the parameters
        potential = ParametricPotential(
            name="mypotential",
            expression="a*x+b",
            parameters={"a": 1.0 * u.g, "b": 1.0 * u.m},
            independent_variables={"x"},
        )
        potential.set_expression(expression="a*x**2+b")

        assert potential.expression == sympy.sympify("a*x**2+b")
        assert potential.parameters == {"a": 1 * u.g, "b": 1 * u.m}

    def test_set_expression_bad(self):
        potential = ParametricPotential(
            name="mypotential",
            expression="a*x+b",
            parameters={"a": 1.0 * u.g, "b": 1.0 * u.m},
            independent_variables={"x"},
        )

        with pytest.raises(ValueError):
            potential.set_expression(expression="a*x**2+b*x+c")

    def test_set_params(self):
        potential = ParametricPotential(
            name="mypotential",
            expression="a*x+b",
            parameters={"a": 1.0 * u.g, "b": 1.0 * u.m},
            independent_variables={"x"},
        )
        potential.set_expression(parameters={"a": 4.0 * u.g, "b": 4.0 * u.m})

        assert potential.expression == sympy.sympify("a*x+b")
        assert potential.parameters == {"a": 4 * u.g, "b": 4 * u.m}

    def test_set_params_bad(self):
        potential = ParametricPotential(
            name="mypotential",
            expression="a*x+b",
            parameters={"a": 1.0 * u.g, "b": 1.0 * u.m},
            independent_variables={"x"},
        )

        with pytest.raises(ValueError):
            potential.set_expression(
                parameters={"aa": 4.0 * u.g, "bb": 4.0 * u.m}
            )

    def test_set_params_partial(self):
        potential = ParametricPotential(
            name="mypotential",
            expression="a*x+b",
            parameters={"a": 1.0 * u.g, "b": 1.0 * u.m},
            independent_variables={"x"},
        )
        potential.set_expression(parameters={"a": 4.0 * u.g})

        assert potential.expression == sympy.sympify("a*x+b")
        assert potential.parameters == {"a": 4 * u.g, "b": 1 * u.m}

    def test_set_expression_and_params(self):
        potential = ParametricPotential(
            name="mypotential",
            expression="a*x+b",
            parameters={"a": 1.0 * u.g, "b": 1.0 * u.m},
            independent_variables={"x"},
        )
        potential.set_expression(
            expression="u*r+v",
            parameters={"u": 1.0 * u.g, "v": 1.0 * u.m},
            independent_variables={"r"},
        )

        assert potential.expression == sympy.sympify("u*r+v")
        assert potential.parameters == {"u": 1 * u.g, "v": 1 * u.m}

    def test_set_expression_and_params_mismatch(self):
        potential = ParametricPotential(
            name="mypotential",
            expression="a*x+b",
            parameters={"a": 1.0 * u.g, "b": 1.0 * u.m},
            independent_variables={"x"},
        )

        with pytest.raises(ValueError):
            potential.set_expression(
                expression="c*x+d",
                parameters={"u": 1.0 * u.g, "v": 1.0 * u.m},
            )

    def test_class_method(self):
        template = PotentialTemplateLibrary()["HarmonicBondPotential"]
        params = {"k": 1.0 * u.kcal / u.nm**2, "r_eq": 1.0 * u.nm}
        harmonic_potential_from_template = ParametricPotential.from_template(
            template, params
        )

        harmonic_potential = ParametricPotential(
            name="HarmonicBondPotential",
            expression="0.5 * k * (r-r_eq)**2",
            independent_variables={"r"},
            parameters=params,
        )

        assert harmonic_potential.name == harmonic_potential_from_template.name
        assert (
            harmonic_potential.expression
            == harmonic_potential_from_template.expression
        )
        assert (
            harmonic_potential.independent_variables
            == harmonic_potential_from_template.independent_variables
        )

    def test_class_method_with_error(self):
        template = object()
        with pytest.raises(TypeError):
            ParametricPotential.from_template(template, parameters=None)

    def test_template_parameterization_dimension_mismatch(self):
        template = PotentialTemplateLibrary()["HarmonicBondPotential"]
        params = {
            "k": 1.0 * u.kcal * u.dimensionless / u.nm,
            "r_eq": 1.0 * u.nm,
        }

        with pytest.raises(AssertionError):
            harmonic_potential_from_template = (
                ParametricPotential.from_template(template, params)
            )

    def test_bondtype_clone(self):
        top = Topology()
        btype = BondType(
            name="ff_255~ff_256",
            expression="a*x+b+c*y",
            independent_variables={"x", "y"},
            parameters={"a": 200 * u.g, "b": 300 * u.K, "c": 400 * u.J},
        )
        btype_clone = btype.clone()

        atom1 = Atom(name="1")
        atom2 = Atom(name="2")
        bond1 = Bond(connection_members=[atom1, atom2])

        atom3 = Atom(name="3")
        atom4 = Atom(name="4")
        bond2 = Bond(connection_members=[atom3, atom4])

        bond1.bond_type = btype
        bond2.bond_type = btype_clone

        top.add_connection(bond1)
        top.add_connection(bond2)
        top.update_topology()

        assert len(top.bond_types) == 2

        btype_dict = btype.model_dump(exclude={"topology", "set_ref"})
        btype_clone_dict = btype_clone.model_dump(
            exclude={"topology", "set_ref"}
        )

        for key, value in btype_dict.items():
            cloned = btype_clone_dict[key]
            assert value == cloned
            if id(value) == id(cloned):
                assert isinstance(value, (str, type(None)))
                assert isinstance(cloned, (str, type(None)))

    def test_sorting(self, parmed_benzene):
        from gmso.external import from_parmed

        top = from_parmed(parmed_benzene)

        labelsList = [
            "bond_types",
            "angle_types",
            "dihedral_types",
            "improper_types",
        ]
        for connection_type in labelsList:
            conn = list(getattr(top, connection_type)())[0]
            assert sort_by_classes(conn) == sort_by_types(conn)
