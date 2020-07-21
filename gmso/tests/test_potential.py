import unyt as u
import sympy
import pytest

from gmso.core.potential import Potential
from gmso.tests.base_test import BaseTest
from gmso.utils.testing import allclose
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.exceptions import GMSOError


class TestPotential(BaseTest):
    def test_new_potential(self):
        new_potential = Potential(
            name='mypotential',
            expression='a*x+b',
            parameters={
                'a': 1.0*u.g,
                'b': 1.0*u.m},
            independent_variables={'x'}
        )
        assert new_potential.name == 'mypotential'
        assert new_potential.expression == sympy.sympify('a*x+b')
        assert allclose(new_potential.parameters['a'], 1.0 * u.g)
        assert allclose(new_potential.parameters['b'], 1.0 * u.m)
        assert new_potential.independent_variables == {sympy.symbols('x')}

    def test_setters(self):
        new_potential = Potential()
        new_potential.name = "SettingName"
        new_potential.set_expression(
            independent_variables=sympy.symbols({'r'}),
            expression='r*sigma*epsilon',
            parameters={
                'sigma': 1 * u.nm,
                'epsilon': 10 * u.Unit('kcal / mol')
            }
        )

        assert new_potential.name == "SettingName"
        assert new_potential.independent_variables == {sympy.symbols('r')}
        assert new_potential.parameters == {
            'sigma': 1 * u.nm,
            'epsilon': 10 * u.Unit('kcal / mol')
        }
        assert new_potential.expression == sympy.sympify('r * sigma * epsilon')

    def test_incorrect_indep_vars(self):
        with pytest.raises(ValueError):
            Potential(expression='x*y', independent_variables='z')

    def test_incorrect_expression(self):
        with pytest.raises(ValueError):
            Potential(name='mytype', expression=4.2)

    def test_expression_consistency(self):
        new_potential = Potential(
            name='mypotential',
            parameters={'x': 1*u.m, 'y': 10*u.m},
            expression='x + y * z',
            independent_variables='z'
        )

        symbol_x, symbol_y, symbol_z = sympy.symbols('x y z')
        correct_expr = sympy.sympify('x+y*z')
        assert new_potential.expression.free_symbols == set([symbol_x, symbol_y, symbol_z])
        assert correct_expr == new_potential.expression

    def test_equivalance(self):
        ref_potential = Potential(
            name='mypotential',
            expression='a*x+b',
            parameters={
                'a': 1.0*u.g,
                'b': 1.0*u.m},
            independent_variables={'x'}
        )
        same_potential = Potential(
            name='mypotential',
            expression='a*x+b',
            parameters={
                'a': 1.0*u.g,
                'b': 1.0*u.m},
            independent_variables={'x'}
        )
        diff_name = Potential(
            name='nopotential',
            expression='a*x+b',
            parameters={
                'a': 1.0*u.g,
                'b': 1.0*u.m},
            independent_variables={'x'}
        )
        diff_expression = Potential(
            name='mypotential',
            expression='a+x*b',
            parameters={
                'a': 1.0*u.g,
                'b': 1.0*u.m},
            independent_variables={'x'}
        )
        diff_params = Potential(
            name='mypotential',
            expression='a*x+b',
            parameters={
                'a': 2.0*u.g,
                'b': 1.0*u.m},
            independent_variables={'x'}
        )

        assert ref_potential    == same_potential
        assert ref_potential    != diff_name
        assert ref_potential    != diff_expression
        assert ref_potential    != diff_params
        assert same_potential   != diff_name
        assert same_potential   != diff_expression
        assert same_potential   != diff_params
        assert diff_name        != diff_expression
        assert diff_expression  != diff_params

    def test_set_expression(self):
        # Try changing the expression but keep the parameters
        potential = Potential(
            name='mypotential',
            expression='a*x+b',
            parameters={
                'a': 1.0*u.g,
                'b': 1.0*u.m},
            independent_variables={'x'}
        )
        potential.set_expression(expression='a*x**2+b')

        assert potential.expression == sympy.sympify('a*x**2+b')
        assert potential.parameters == {'a': 1*u.g, 'b': 1*u.m}

    def test_set_expression_bad(self):
        potential = Potential(
            name='mypotential',
            expression='a*x+b',
            parameters={
                'a': 1.0*u.g,
                'b': 1.0*u.m},
            independent_variables={'x'}
        )

        with pytest.raises(ValueError):
            potential.set_expression(expression='a*x**2+b*x+c')

    def test_set_params(self):
        potential = Potential(
            name='mypotential',
            expression='a*x+b',
            parameters={
                'a': 1.0*u.g,
                'b': 1.0*u.m},
            independent_variables={'x'}
        )
        potential.set_expression(parameters={'a': 4.0*u.g, 'b': 4.0*u.m})

        assert potential.expression == sympy.sympify('a*x+b')
        assert potential.parameters == {'a': 4*u.g, 'b': 4*u.m}

    def test_set_params_bad(self):
        potential = Potential(
            name='mypotential',
            expression='a*x+b',
            parameters={
                'a': 1.0*u.g,
                'b': 1.0*u.m},
            independent_variables={'x'}
        )

        with pytest.raises(ValueError):
            potential.set_expression(parameters={'aa': 4.0*u.g, 'bb': 4.0*u.m})

    def test_set_params_partial(self):
        potential = Potential(
            name='mypotential',
            expression='a*x+b',
            parameters={
                'a': 1.0*u.g,
                'b': 1.0*u.m},
            independent_variables={'x'}
        )
        potential.set_expression(parameters={'a': 4.0*u.g})

        assert potential.expression == sympy.sympify('a*x+b')
        assert potential.parameters == {'a': 4*u.g, 'b': 1*u.m}

    def test_set_expression_and_params(self):
        potential = Potential(
            name='mypotential',
            expression='a*x+b',
            parameters={
                'a': 1.0*u.g,
                'b': 1.0*u.m},
            independent_variables={'x'}
        )
        potential.set_expression(
            expression='u*r+v',
            parameters={
                'u': 1.0*u.g,
                'v': 1.0*u.m},
            independent_variables={'r'}
        )

        assert potential.expression == sympy.sympify('u*r+v')
        assert potential.parameters == {'u': 1*u.g, 'v': 1*u.m}

    def test_set_expression_and_params_mismatch(self):
        potential = Potential(
            name='mypotential',
            expression='a*x+b',
            parameters={
                'a': 1.0*u.g,
                'b': 1.0*u.m},
            independent_variables={'x'}
        )

        with pytest.raises(ValueError):
            potential.set_expression(
                expression='c*x+d',
                parameters={
                    'u': 1.0*u.g,
                    'v': 1.0*u.m},
            )

    def test_class_method(self):
        template = PotentialTemplateLibrary()['HarmonicBondPotential']
        params = {'k': 1.0 * u.dimensionless,
                  'r_eq': 1.0 * u.dimensionless}
        harmonic_potential_from_template = Potential.from_template(template, params)
        harmonic_potential = Potential(name='HarmonicBondPotential',
                                       expression='0.5 * k * (r-r_eq)**2',
                                       independent_variables={'r'},
                                       parameters=params,
                                       template=False)
        assert harmonic_potential.name == harmonic_potential_from_template.name
        assert harmonic_potential.expression == harmonic_potential_from_template.expression
        assert harmonic_potential.independent_variables == harmonic_potential_from_template.independent_variables

    def test_class_method_with_error(self):
        template = object()
        with pytest.raises(GMSOError):
            Potential.from_template(template, parameters=None)
