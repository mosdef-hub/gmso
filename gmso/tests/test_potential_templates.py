import os
import glob
import pytest

import sympy

from gmso.lib.potential_templates import PotentialTemplates, JSON_DIR
from gmso.exceptions import GMSOError
from gmso.tests.base_test import BaseTest


class TestPotentialTemplates(BaseTest):
    @pytest.fixture
    def templates(self):
        return PotentialTemplates()

    def test_singleton_behavior(self, templates):
        assert id(templates) == id(PotentialTemplates())
        assert id(PotentialTemplates()) == id(PotentialTemplates())

    def test_lennard_jones_potential(self, templates):
        lennard_jones_potential = templates['LennardJonesPotential']
        assert lennard_jones_potential.name == 'LennardJonesPotential'
        assert lennard_jones_potential.expression == sympy.sympify('4*epsilon*((sigma/r)**12 - (sigma/r)**6)')
        assert lennard_jones_potential.independent_variables == {sympy.sympify('r')}

    def test_mie_potential(self, templates):
        mie_potential = templates['MiePotential']
        assert mie_potential.name == 'MiePotential'
        assert mie_potential.expression == sympy.sympify(
            '(n/(n-m)) * (n/m)**(m/(n-m)) * epsilon * ((sigma/r)**n - (sigma/r)**m)')
        assert mie_potential.independent_variables == {sympy.sympify('r')}

    def test_opls_torsion_potential(self, templates):
        opls_torsion_potential = templates['OPLSTorsionPotential']
        assert opls_torsion_potential.name == 'OPLSTorsionPotential'
        assert opls_torsion_potential.expression == sympy.sympify('0.5 * k0 + 0.5 * k1 * (1 + cos(phi)) + '
                                                                  '0.5 * k2 * (1 - cos(2*phi)) +'
                                                                  ' 0.5 * k3 * (1 + cos(3*phi))'
                                                                  ' + 0.5 * k4 * (1 - cos(4*phi))')
        assert opls_torsion_potential.independent_variables == {sympy.sympify('phi')}

    def test_periodic_torsion_potential(self, templates):
        periodic_torsion_potential = templates['PeriodicTorsionPotential']
        assert periodic_torsion_potential.name == 'PeriodicTorsionPotential'
        assert periodic_torsion_potential.expression == sympy.sympify('k * (1 + cos(n * phi - phi_eq))**2')
        assert periodic_torsion_potential.independent_variables == {sympy.sympify('phi')}

    def test_ryckaert_bellemans_torsion_potential(self, templates):
        ryckaert_bellemans_torsion_potential = templates['RyckaertBellemansTorsionPotential']
        assert ryckaert_bellemans_torsion_potential.name == 'RyckaertBellemansTorsionPotential'
        assert ryckaert_bellemans_torsion_potential.expression == sympy.sympify('c0 * cos(phi)**0 + c1 * cos(phi)**1 +'
                                                                                ' c2 * cos(phi)**2 + c3 * cos(phi)**3 +'
                                                                                ' c4 * cos(phi)**4 + c5 * cos(phi)**5')
        assert ryckaert_bellemans_torsion_potential.independent_variables == {sympy.sympify('phi')}

    def test_harmonic_torsion_potential(self, templates):
        harmonic_torsion_potential = templates['HarmonicTorsionPotential']
        assert harmonic_torsion_potential.name == 'HarmonicTorsionPotential'
        assert harmonic_torsion_potential.expression == sympy.sympify('0.5 * k * (phi - phi_eq)**2')
        assert harmonic_torsion_potential.independent_variables == {sympy.sympify('phi')}

    def test_harmonic_improper_potential(self, templates):
        harmonic_improper_potential = templates['HarmonicImproperPotential']
        assert harmonic_improper_potential.name == 'HarmonicImproperPotential'
        assert harmonic_improper_potential.expression == sympy.sympify('0.5 * k * (phi - phi_eq)**2')
        assert harmonic_improper_potential.independent_variables == {sympy.sympify('phi')}

    def test_harmonic_bond_potential(self, templates):
        harmonic_bond_potential = templates['HarmonicBondPotential']
        assert harmonic_bond_potential.name == 'HarmonicBondPotential'
        assert harmonic_bond_potential.expression == sympy.sympify('0.5 * k * (r-r_eq)**2')
        assert harmonic_bond_potential.independent_variables == {sympy.sympify('r')}

    def test_harmonic_angle_potential(self, templates):
        harmonic_bond_potential = templates['HarmonicAnglePotential']
        assert harmonic_bond_potential.name == 'HarmonicAnglePotential'
        assert harmonic_bond_potential.expression == sympy.sympify('0.5 * k * (theta-theta_eq)**2')
        assert harmonic_bond_potential.independent_variables == {sympy.sympify('theta')}

    def test_buckingham_potential(self, templates):
        buckingham_potential = templates['BuckinghamPotential']
        assert buckingham_potential.name == 'BuckinghamPotential'
        assert buckingham_potential.expression == sympy.sympify('a*exp(-b*r) - c*r**-6')
        assert buckingham_potential.independent_variables == sympy.sympify({'r'})

    def test_save_templates(self, templates):
        PotentialTemplates.save_potential_template('TestTemplate',
                                                   {
                                                       'name': 'TestTemplate',
                                                       'expression': 'a+b+c',
                                                       'independent_variables': 'c,a'
                                                   },
                                                   update=True,
                                                   user_template=False
                                                   )

        assert os.path.exists(os.path.join(JSON_DIR, 'TestTemplate.json'))
        assert templates['TestTemplate']

    def test_available_template(self, templates):
        names = templates.get_available_template_names()
        assert isinstance(names, tuple)
        assert len(names) == len(glob.glob(os.path.join(JSON_DIR, '*.json')))

    def test_name_mismatch_in_templates(self):
        wrong_ref_dict = {
                'name': 'NameMismatchTemplate',
                'expression': 'a+b+c',
                'independent_variables': 'c,a'
            }
        with pytest.raises(GMSOError):
            PotentialTemplates.save_potential_template('NameMismatch_Template', wrong_ref_dict, update=False)

    def test_user_defined_potential(self, templates):
        user_potential_dict = {
            'name': 'UserTemplate',
            'expression': 'a*b*c',
            'independent_variables': 'b'
        }
        PotentialTemplates.save_potential_template('UserTemplate', user_potential_dict, update=True, user_template=True)
        assert os.path.exists(os.path.join(os.path.expanduser('~'),
                                           '.gmso',
                                           'potential_templates',
                                           'UserTemplate.json'))
        templates.load_user_templates()
        assert templates['UserTemplate'].name == 'UserTemplate'
        assert templates['UserTemplate'].expression == sympy.sympify('a*b*c')
        assert templates['UserTemplate'].independent_variables == {sympy.sympify('b')}
