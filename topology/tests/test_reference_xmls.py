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

    def test_spce_xml(self):
        spce = ForceField(get_path('spce.xml'))

        assert len(spce.atom_types) == 2
        assert len(spce.bond_types) == 1
        assert len(spce.angle_types) == 1
        assert len(spce.dihedral_types) == 0

        # Store expected expressions in list
        ref_exprs = [sympy.sympify(expr) for expr in [
            "4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            "0.5 * k * (r-r_eq)**2",
            "0.5 * k * (theta-theta_eq)**2",
            ]
        ]

        assert allclose(spce.atom_types['opls_116'].charge, -0.8476 * u.elementary_charge)
        assert allclose(spce.atom_types['opls_117'].charge, 0.4238 * u.elementary_charge)

        assert sympy.simplify(spce.atom_types['opls_116'].expression - ref_exprs[0]) == 0
        assert sympy.simplify(spce.atom_types['opls_117'].expression - ref_exprs[0]) == 0
        assert sympy.simplify(spce.bond_types['opls_116~opls_117'].expression - ref_exprs[1]) == 0
        assert sympy.simplify(spce.angle_types['opls_117~opls_116~opls_117'].expression - ref_exprs[2]) == 0

        assert allclose(spce.atom_types['opls_116'].parameters['sigma'], 0.316557 * u.nm)
        assert allclose(spce.atom_types['opls_116'].parameters['epsilon'], 0.650194 * u.Unit('kJ/mol'))
        assert allclose(spce.atom_types['opls_117'].parameters['sigma'], 0.1 * u.nm)
        assert allclose(spce.atom_types['opls_117'].parameters['epsilon'], 0.0 * u.Unit('kJ/mol'))

        assert allclose(spce.bond_types['opls_116~opls_117'].parameters['r_eq'], 0.1 * u.nm)
        assert allclose(spce.bond_types['opls_116~opls_117'].parameters['k'], 345000.0 * u.Unit('kJ/(mol*nm**2)'))

        assert allclose(spce.angle_types['opls_117~opls_116~opls_117'].parameters['theta_eq'], 109.47 * u.degree)
