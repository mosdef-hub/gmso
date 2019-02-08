import numpy as np
import sympy
import pytest

from topology.core.atom_type import AtomType


def test_new_atom_type():
    new_type = AtomType(name='mytype', charge=1.0, 
            parameters={'sigma':1, 'epsilon':10})
    assert new_type.name == 'mytype'
    assert np.isclose(new_type.charge, 1.0)
    assert np.isclose(new_type.parameters['sigma'], 1)
    assert np.isclose(new_type.parameters['epsilon'], 10)

def test_setters():
    new_type = AtomType()
    new_type.name = "SettingName"
    new_type.charge = -1.0
    new_type.parameters = {'sigma':100, 'epsilon':400}
    assert new_type.name == "SettingName"
    assert np.isclose(new_type.charge , -1.0)
    assert new_type.parameters == {'sigma':100, 'epsilon':400}

def test_incorrect_nb_function():
    with pytest.raises(ValueError):
        new_type = AtomType(name='mytype', charge=1.0, nb_function=4.2)

def test_nb_function_consistency():
    new_type = AtomType(name='mytype', charge=1.0, 
            parameters={'x':1, 'y':10}, nb_function='x + y')
    symbol_x, symbol_y = sympy.symbols('x y')
    assert new_type.nb_function.free_symbols == set([symbol_x, symbol_y])

def test_equivalance():
    first_type = AtomType(name='mytype', charge=1.0, 
            parameters={'sigma':1, 'epsilon':10})
    second_type = AtomType(name='mytype', charge=1.0, 
            parameters={'sigma':1, 'epsilon':10})
    different_type = AtomType(name='difftype', charge=4.0, 
            parameters={'sigma':1, 'epsilon':10})

    
    assert first_type == second_type
    assert first_type != different_type
