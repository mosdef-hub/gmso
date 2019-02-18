import unyt as u
import numpy as np
import sympy
import pytest

from topology.core.atom_type import AtomType
from topology.tests.base_test import BaseTest
from topology.testing.utils import allclose


class TestAtomType(BaseTest):

    def test_new_atom_type(self, charge):
        new_type = AtomType(name='mytype', charge=charge,
                parameters={'sigma':1 * u.nm, 
                    'epsilon':10 * u.Unit('kcal / mol')})
        assert new_type.name == 'mytype'
        assert allclose(new_type.charge, charge)
        assert allclose(new_type.parameters['sigma'], 1 * u.nm)
        assert allclose(new_type.parameters['epsilon'], 10 * u.Unit('kcal / mol'))
    
    def test_setters(self, charge):
        new_type = AtomType(self)
        new_type.name = "SettingName"
        new_type.charge = -1.0 * charge
        assert new_type.name == "SettingName"
        assert allclose(new_type.charge , -1.0 * charge)
    
    def test_incorrect_nb_function(self, charge):
        with pytest.raises(ValueError):
            new_type = AtomType(name='mytype', charge=charge, nb_function=4.2)
    
    def test_nb_function_consistency(self, charge):
        # Test nb-func symbol consistency with parameter consistency in init
        new_type = AtomType(name='mytype', charge=charge,
                parameters={'x':1, 'y':10}, nb_function='x + y')
        symbol_x, symbol_y = sympy.symbols('x y')
        correct_expr = sympy.sympify('x+y')
        assert new_type.nb_function.free_symbols == set([symbol_x, symbol_y])
        assert correct_expr == new_type.nb_function
    
    def test_equivalance(self, charge):
        first_type = AtomType(name='mytype', charge=charge,
                parameters={'sigma':1, 'epsilon':10})
        same_type = AtomType(name='mytype', charge=charge,
                parameters={'sigma':1, 'epsilon':10})
        different_name = AtomType(name='difftype', charge=charge,
                parameters={'sigma':1, 'epsilon':10})
        different_charge = AtomType(name='mytype', charge=4.0 * charge, 
                parameters={'sigma':1, 'epsilon':10})
        different_function = AtomType(name='mytype', charge=charge,
                parameters={'sigma':1, 'epsilon':10},
                nb_function='sigma * epsilon')
        different_params = AtomType(name='mytype', charge=charge,
                parameters={'sigma':42, 'epsilon':100000})
                
    
        assert first_type == same_type
        assert first_type != different_name
        assert first_type != different_charge
        assert first_type != different_function
        assert first_type != different_params
    
    def test_set_nb_func(self, charge):
        # Try changing the nonbonded function, but keep the parameters
        first_type = AtomType(name='mytype', charge=charge,
                parameters={'sigma':1, 'epsilon':10})
        first_type.set_nb_function(function='sigma + epsilon')
        correct_expr = sympy.sympify('sigma + epsilon')
        assert first_type.nb_function == correct_expr
        assert first_type.parameters == {'sigma':1, 'epsilon':10}
    
    def test_set_nb_func_bad(self):
        # Try changing the nonbonded function, keep the parameters, 
        # but the nb function uses different symbols
        first_type = AtomType(nb_function='sigma*epsilon', 
                parameters={'sigma':1, 'epsilon':100})
        with pytest.raises(ValueError):
            first_type.set_nb_function(function='a + b')
    
    def test_set_nb_func_params(self, charge):
        # Try changing the nb parameters, but keeping the function
        first_type = AtomType(name='mytype', charge=charge,
                nb_function='sigma * epsilon',
                parameters={'sigma':1, 'epsilon':10})
        first_type.set_nb_function(parameters={'sigma':42, 'epsilon':24})
        correct_expr = sympy.sympify('sigma * epsilon')
        assert first_type.nb_function == correct_expr
        assert first_type.parameters == {'sigma':42, 'epsilon':24}
    
    def test_set_nb_func_params_bad(self):
        # Try changing the parameters, keep the function,
        # but the new parameters use different symbols
        first_type = AtomType(nb_function='sigma*epsilon',
                parameters={'sigma':1, 'epsilon':10})
        with pytest.warns(UserWarning):
            first_type.set_nb_function(parameters={'a':1, 'b':10})
    
    def test_set_nb_func_params_partial(self):
        # Try changing the parameters, keep the function,
        # but change only one symbol
        first_type = AtomType(nb_function='sigma*epsilon',
                parameters={'sigma':1, 'epsilon':10})
        first_type.set_nb_function(parameters={'sigma':42})
        correct_expr = sympy.sympify('sigma*epsilon')
        assert first_type.parameters == {'sigma':42, 'epsilon':10}
        assert first_type.nb_function == correct_expr
    
    def test_set_nb_func_params_both_correct(self):
        # Try correctly changing both the nb function and parameters
        first_type = AtomType(nb_function='sigma*epsilon',
                parameters={'sigma':1, 'epsilon':10})
        first_type.set_nb_function(function='a+b', parameters={'a':100, 'b':42})
        correct_expr = sympy.sympify('a+b')
        correct_params = {'a':100, 'b':42}
        assert first_type.nb_function == correct_expr
        assert first_type.parameters == correct_params
    
    def test_set_nb_func_params_both_incorrect(self):
        # Try incorrectly changing both the nb function and the parameters
        first_type = AtomType(nb_function='sigma*epsilon',
                parameters={'sigma':1, 'epsilon':10})
        with pytest.raises(ValueError):
            first_type.set_nb_function(function='a+b', parameters={'c':100, 'd':42})
