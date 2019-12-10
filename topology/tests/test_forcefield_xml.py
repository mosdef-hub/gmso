import pytest
from sympy import sympify
import unyt as u

from topology.forcefield import ForceField
from topology.tests.utils import get_path
from topology.tests.base_test import BaseTest


class TestForceFieldFromXML(BaseTest):

    @pytest.fixture
    def ff(self):
        return ForceField.from_xml(get_path('ff-example0.xml'))

    def test_ff_name_version_from_xml(self, ff):
        assert ff.name == 'ForceFieldOne'
        assert ff.version == '0.4.1'

    def test_ff_atomtypes_from_xml(self, ff):
        assert len(ff.atom_types) == 2
        assert 'Ar' in ff.atom_types
        assert 'Xe' in ff.atom_types

        assert sympify('r') in ff.atom_types['Ar'].independent_variables
        assert ff.atom_types['Ar'].parameters['A'] == u.unyt_quantity(0.1, u.kcal / u.mol)
        assert ff.atom_types['Ar'].parameters['B'] == u.unyt_quantity(4.0, u.nm)
        assert ff.atom_types['Ar'].parameters['C'] == u.unyt_quantity(0.5, u.kcal / u.mol * u.nm ** 6)
        assert ff.atom_types['Ar'].mass == u.unyt_quantity(39.948, u.amu)
        assert ff.atom_types['Ar'].charge == u.unyt_quantity(0.0, u.coulomb)
        assert ff.atom_types['Ar'].description == 'Argon atom'
        assert ff.atom_types['Ar'].definition == 'Ar'
        assert ff.atom_types['Ar'].expression == sympify('(A*exp(-B/r) - C/r**6)')

        assert sympify('r') in ff.atom_types['Xe'].independent_variables
        assert 'A' in ff.atom_types['Xe'].parameters
        assert ff.atom_types['Xe'].parameters['A'] == u.unyt_quantity(0.2, u.kcal / u.mol)
        assert ff.atom_types['Xe'].parameters['B'] == u.unyt_quantity(5.0, u.nm)
        assert ff.atom_types['Xe'].parameters['C'] == u.unyt_quantity(0.3, u.kcal / u.mol * u.nm ** 6)
        assert ff.atom_types['Xe'].mass == u.unyt_quantity(131.293, u.amu)
        assert ff.atom_types['Xe'].charge == u.unyt_quantity(0.0, u.coulomb)
        assert ff.atom_types['Xe'].description == 'Xenon atom'
        assert ff.atom_types['Xe'].definition == 'Xe'
        assert ff.atom_types['Xe'].expression == sympify('(A*exp(-B/r) - C/r**6)')

    def test_ff_bondtypes_from_xml(self, ff):
        assert len(ff.bond_types) == 2
        assert 'BondType1' in ff.bond_types
        assert 'BondType2' in ff.bond_types

        assert sympify('r') in ff.bond_types['BondType1'].independent_variables
        assert ff.bond_types['BondType1'].parameters['r_eq'] == u.unyt_quantity(10.0, u.nm)
        assert ff.bond_types['BondType1'].parameters['k'] == u.unyt_quantity(10000, u.kJ / u.mol)
        assert ff.bond_types['BondType1'].member_types == ['Ar', 'Ar']

        assert sympify('r') in ff.bond_types['BondType2'].independent_variables
        assert ff.bond_types['BondType2'].parameters['r_eq'] == u.unyt_quantity(10.0, u.nm)
        assert ff.bond_types['BondType2'].parameters['k'] == u.unyt_quantity(20000, u.kJ / u.mol)
        assert ff.bond_types['BondType2'].member_types == ['Xe', 'Xe']

    def test_ff_angletypes_from_xml(self, ff):
        assert len(ff.angle_types) == 2
        assert 'AngleType1' in ff.angle_types
        assert 'AngleType2' in ff.angle_types

        assert sympify('r') in ff.angle_types['AngleType1'].independent_variables
        assert ff.angle_types['AngleType1'].parameters['r_eq'] == u.unyt_quantity(10.0, u.nm)
        assert ff.angle_types['AngleType1'].parameters['z'] == u.unyt_quantity(100, u.kJ / u.mol)
        assert ff.angle_types['AngleType1'].member_types == ['Ar', 'Ar', 'Ar']

        assert sympify('r') in ff.angle_types['AngleType2'].independent_variables
        assert ff.angle_types['AngleType2'].parameters['r_eq'] == u.unyt_quantity(10.0, u.nm)
        assert ff.angle_types['AngleType2'].parameters['z'] == u.unyt_quantity(20, u.kJ / u.mol)
        assert ff.angle_types['AngleType2'].member_types == ['Xe', 'Xe', 'Xe']

    def test_ff_dihedraltypes_from_xml(self, ff):
        assert len(ff.dihedral_types) == 2
        assert 'DihedralType1' in ff.dihedral_types
        assert 'DihedralType2' in ff.dihedral_types

        assert sympify('r') in ff.dihedral_types['DihedralType1'].independent_variables
        assert ff.dihedral_types['DihedralType1'].parameters['r_eq'] == u.unyt_quantity(10.0, u.nm)
        assert ff.dihedral_types['DihedralType1'].parameters['z'] == u.unyt_quantity(100, u.kJ / u.mol)
        assert ff.dihedral_types['DihedralType1'].member_types == ['Ar', 'Ar', 'Ar', 'Ar']

        assert sympify('r') in ff.dihedral_types['DihedralType2'].independent_variables
        assert ff.dihedral_types['DihedralType2'].parameters['r_eq'] == u.unyt_quantity(10.0, u.nm)
        assert ff.dihedral_types['DihedralType2'].parameters['z'] == u.unyt_quantity(20, u.kJ / u.mol)
        assert ff.dihedral_types['DihedralType2'].member_types == ['Xe', 'Xe', 'Xe', 'Xe']
