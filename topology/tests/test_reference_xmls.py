import numpy as np
import sympy
import unyt as u

from topology.forcefield import ForceField
from topology.utils.testing import allclose
from topology.tests.utils import get_path
from topology.tests.base_test import BaseTest


class TestForceFieldFromXML(BaseTest):
    def test_carbon_force_field(self):
        carbon = ForceField(get_path('carbon.xml'))

        assert len(carbon.atom_types) == 1
        assert len(carbon.bond_types) == 1
        assert len(carbon.angle_types) == 1
        assert len(carbon.dihedral_types) == 1

        # Store expected expressions in list
        ref_exprs = [sympy.sympify(expr) for expr in [
            "4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            "0.5 * k * (r-r_eq)**2",
            "0.5 * k * (theta-theta_eq)**2",
            "k * (1 + cos(n * theta - theta_0))",
            ]
        ]

        assert carbon.atom_types['C'].charge.value == 0

        assert sympy.simplify(carbon.atom_types['C'].expression - ref_exprs[0]) == 0
        assert sympy.simplify(carbon.bond_types['C~C'].expression - ref_exprs[1]) == 0
        assert sympy.simplify(carbon.angle_types['C~C~C'].expression - ref_exprs[2]) == 0
        assert sympy.simplify(carbon.dihedral_types['*~C~C~*'].expression - ref_exprs[3]) == 0

        assert allclose(carbon.atom_types['C'].parameters['sigma'], 0.339966950842 * u.nm)
        assert allclose(carbon.atom_types['C'].parameters['epsilon'], 0.359824 * u.Unit('kJ/mol'))

        assert allclose(carbon.bond_types['C~C'].parameters['r_eq'], 0.1324 * u.nm)
        assert allclose(carbon.bond_types['C~C'].parameters['k'], 493460.96 * u.Unit('kJ/(mol*nm**2)'))

        assert allclose(carbon.angle_types['C~C~C'].parameters['theta_eq'], 2.12598556185 * u.radian)
        assert allclose(carbon.angle_types['C~C~C'].parameters['k'], 584.42112 * u.Unit('kJ/(mol*rad**2)'))

        assert allclose(carbon.dihedral_types['*~C~C~*'].parameters['k'], 27.8236 * u.Unit('kJ/mol'))
        assert allclose(carbon.dihedral_types['*~C~C~*'].parameters['n'], 2 * u.dimensionless)
        assert allclose(carbon.dihedral_types['*~C~C~*'].parameters['theta_0'], np.pi * u.radian)

    def test_tip3p_force_field(self):
        water = ForceField(get_path('tip3p.xml'))
        assert len(water.atom_types) == 2
        assert len(water.bond_types) == 1
        assert len(water.angle_types) == 1
        assert len(water.dihedral_types) == 0

        # Store expected expressions in list
        ref_exprs = [sympy.sympify(expr) for expr in [
            "4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            "0.5 * k * (r-r_eq)**2",
            "0.5 * k * (theta-theta_eq)**2",
            ]
        ]

        assert water.atom_types['opls_111'].charge.value == -0.834
        assert water.atom_types['opls_112'].charge.value == 0.417

        assert sympy.simplify(water.atom_types['opls_111'].expression - ref_exprs[0]) == 0
        assert sympy.simplify(water.atom_types['opls_112'].expression - ref_exprs[0]) == 0
        assert sympy.simplify(water.bond_types['opls_111~opls_112'].expression - ref_exprs[1]) == 0
        assert sympy.simplify(water.angle_types['opls_112~opls_111~opls_112'].expression - ref_exprs[2]) == 0

        assert allclose(water.atom_types['opls_111'].parameters['sigma'], 0.315061 * u.nm)
        assert allclose(water.atom_types['opls_111'].parameters['epsilon'], 0.636386 * u.Unit('kJ/mol'))
        assert allclose(water.atom_types['opls_112'].parameters['sigma'], 1.0 * u.nm)
        assert allclose(water.atom_types['opls_112'].parameters['epsilon'], 0.0 * u.Unit('kJ/mol'))

        assert allclose(water.bond_types['opls_111~opls_112'].parameters['r_eq'], 0.09572 * u.nm)
        assert allclose(water.bond_types['opls_111~opls_112'].parameters['k'], 502416.0 * u.Unit('kJ/(mol*nm**2)'))

        assert allclose(water.angle_types['opls_112~opls_111~opls_112'].parameters['theta_eq'], 1.824218134 * u.radian)
        assert allclose(water.angle_types['opls_112~opls_111~opls_112'].parameters['k'], 682.02 * u.Unit('kJ/(mol*rad**2)'))
