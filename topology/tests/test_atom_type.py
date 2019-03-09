import unyt as u
import numpy as np
import sympy
import pytest

from topology.core.atom_type import AtomType
from topology.tests.base_test import BaseTest
from topology.testing.utils import allclose


class TestAtomType(BaseTest):
    def test_new_atom_type(self, charge, mass):
        new_type = AtomType(name='mytype', charge=charge, mass=mass,
                parameters={'sigma': 1 * u.nm,
                    'epsilon': 10 * u.Unit('kcal / mol')},
                independent_variables={'r'})
        assert new_type.name == 'mytype'
        assert allclose(new_type.charge, charge)
        assert allclose(new_type.parameters['sigma'], 1 * u.nm)
        assert allclose(new_type.parameters['epsilon'], 10 * u.Unit('kcal / mol'))
        assert allclose(new_type.mass, mass)

    def test_setters(self, charge, mass):
        new_type = AtomType(self)
        new_type.name = "SettingName"
        new_type.charge = -1.0 * charge
        new_type.mass = 1 * mass
        new_type.independent_variables = 'r'
        new_type.parameters = {'sigma': 1 * u.nm,
                               'epsilon': 10 * u.Unit('kcal / mol')}
        new_type.nb_function = 'r * sigma * epsilon'
        assert new_type.name == "SettingName"
        assert allclose(new_type.charge, -1.0 * charge)
        assert allclose(new_type.mass, 1 * mass)
        assert new_type.independent_variables == {sympy.symbols('r')}
        assert new_type.parameters == {'sigma': 1 * u.nm,
                                      'epsilon': 10 * u.Unit('kcal / mol')}
        assert new_type.nb_function == sympy.sympify('r * sigma * epsilon')

    def test_incorrect_nb_function(self, charge):
        with pytest.raises(ValueError):
            new_type = AtomType(name='mytype', charge=charge, nb_function=4.2)

    def test_nb_function_consistency(self, charge):
        # Test nb-func symbol consistency with parameter consistency in init
        new_type = AtomType(name='mytype', charge=charge,
                parameters={'x': 1*u.m, 'y': 10*u.m}, nb_function='x + y')
        symbol_x, symbol_y = sympy.symbols('x y')
        correct_expr = sympy.sympify('x+y')
        assert new_type.nb_function.free_symbols == set([symbol_x, symbol_y])
        assert correct_expr == new_type.nb_function

    def test_equivalance(self, charge):
        first_type = AtomType(name='mytype', charge=charge,
                parameters={'sigma': 1*u.m, 'epsilon': 10*u.m})
        same_type = AtomType(name='mytype', charge=charge,
                parameters={'sigma': 1*u.m, 'epsilon': 10*u.m})
        different_name = AtomType(name='difftype', charge=charge,
                parameters={'sigma': 1*u.m, 'epsilon': 10*u.m})
        different_charge = AtomType(name='mytype', charge=4.0 * charge,
                parameters={'sigma': 1*u.m, 'epsilon': 10*u.m})
        different_function = AtomType(name='mytype', charge=charge,
                parameters={'sigma': 1*u.m, 'epsilon': 10*u.m},
                nb_function='sigma * epsilon')
        different_params = AtomType(name='mytype', charge=charge,
                parameters={'sigma': 42*u.m, 'epsilon': 100000*u.m})
        different_mass = AtomType(name='mytype', charge=charge, mass=5*u.kg/u.mol,
                parameters={'sigma': 1*u.m, 'epsilon': 10*u.m})

        assert first_type == same_type
        assert first_type != different_name
        assert first_type != different_charge
        assert first_type != different_function
        assert first_type != different_params
        assert first_type != different_mass

    def test_set_nb_func(self, charge):
        # Try changing the nonbonded function, but keep the parameters
        first_type = AtomType(name='mytype', charge=charge,
                parameters={'sigma': 1*u.m, 'epsilon': 10*u.m},
                independent_variables='r')
        first_type.set_nb_function(function='r * (sigma + epsilon)')
        correct_expr = sympy.sympify('r * (sigma + epsilon)')
        assert first_type.nb_function == correct_expr
        assert first_type.parameters == {'sigma': 1*u.m, 'epsilon': 10*u.m}

    def test_set_nb_func_bad(self):
        # Try changing the nonbonded function, keep the parameters,
        # but the nb function uses different symbols
        first_type = AtomType(nb_function='sigma*epsilon',
                parameters={'sigma': 1*u.m, 'epsilon': 100*u.m})
        with pytest.raises(ValueError):
            first_type.set_nb_function(function='a + b')

    def test_set_nb_func_params(self, charge):
        # Try changing the nb parameters, but keeping the function
        first_type = AtomType(name='mytype', charge=charge,
                nb_function='r * sigma * epsilon',
                parameters={'sigma': 1*u.m, 'epsilon': 10*u.m},
                independent_variables='r')
        first_type.set_nb_function(parameters={'sigma': 42*u.m, 'epsilon': 24*u.m})
        correct_expr = sympy.sympify('r * sigma * epsilon')
        assert first_type.nb_function == correct_expr
        assert first_type.parameters == {'sigma': 42, 'epsilon': 24}

    def test_set_nb_func_params_bad(self):
        # Try changing the parameters, keep the function,
        # but the new parameters use different symbols
        first_type = AtomType(nb_function='r*sigma*epsilon',
                parameters={'sigma': 1*u.m, 'epsilon': 10*u.J},
                independent_variables={'r'})

        with pytest.raises(ValueError):
            first_type.set_nb_function(parameters={'a': 1*u.g, 'b': 10*u.m})

    def test_set_nb_func_params_partial(self):
        # Try changing the parameters, keep the function,
        # but change only one symbol
        first_type = AtomType(nb_function='r*sigma*epsilon',
                parameters={'sigma': 1*u.m, 'epsilon': 10*u.J},
                independent_variables={'r'})
        first_type.set_nb_function(parameters={'sigma': 42*u.m})
        correct_expr = sympy.sympify('r*sigma*epsilon')
        assert first_type.parameters == {'sigma': 42*u.m, 'epsilon': 10*u.J}
        assert first_type.nb_function == correct_expr

    def test_set_nb_func_params_both_correct(self):
        # Try correctly changing both the nb function and parameters
        first_type = AtomType(nb_function='sigma*epsilon',
                parameters={'sigma': 1*u.m, 'epsilon': 10*u.J})
        first_type.set_nb_function(function='a+b*x', parameters={'a': 100*u.m, 'b': 42*u.J}, independent_variables='x')
        correct_expr = sympy.sympify('a+b*x')
        correct_params = {'a': 100*u.m, 'b': 42*u.J}
        assert first_type.nb_function == correct_expr
        assert first_type.parameters == correct_params

    def test_set_nb_func_params_both_incorrect(self):
        # Try incorrectly changing both the nb function and the parameters
        first_type = AtomType(nb_function='r*sigma*epsilon',
                parameters={'sigma': 1*u.m, 'epsilon': 10*u.J},
                independent_variables='r')
        with pytest.raises(ValueError):
            first_type.set_nb_function(function='a*x+b',
                parameters={'c': 100*u.year, 'd': 42*u.newton},
                independent_variables='x')
