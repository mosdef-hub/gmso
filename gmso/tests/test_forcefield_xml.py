import lxml
import pytest
from sympy import sympify
import unyt as u

from lxml.etree import DocumentInvalid

from gmso.core.forcefield import ForceField
from gmso.tests.utils import get_path
from gmso.tests.base_test import BaseTest
from gmso.exceptions import ForceFieldParseError


class TestForceFieldFromXML(BaseTest):

    @pytest.fixture
    def ff(self):
        return ForceField(get_path('ff-example0.xml'))

    @pytest.fixture
    def named_groups_ff(self):
        return ForceField(get_path('ff-example1.xml'))

    def test_ff_name_version_from_xml(self, ff):
        assert ff.name == 'ForceFieldOne'
        assert ff.version == '0.4.1'

    def test_scaling_factors_from_xml(self, ff):
        assert ff.scaling_factors['nonBonded14Scale'] == 0.67
        assert ff.scaling_factors['electrostatics14Scale'] == 0.5

    @pytest.mark.parametrize('unit_name,unit_value', [('energy', u.Unit(u.K*u.kb)),
                                                      ('mass', u.gram/u.mol), ('temperature', u.K),
                                                      ('charge', u.coulomb), ('angle', u.rad),
                                                      ('time', u.ps), ('distance', u.nm)])
    def test_units_from_xml(self, ff, unit_name, unit_value):
        assert len(ff.units.keys()) == 7
        assert ff.units[unit_name] == unit_value

    def test_ff_atomtypes_from_xml(self, ff):
        assert len(ff.atom_types) == 3
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

        assert ff.atom_types['Li'].charge == u.unyt_quantity(1.0, u.coulomb)

    def test_ff_bondtypes_from_xml(self, ff):
        assert len(ff.bond_types) == 2
        assert 'Ar~Ar' in ff.bond_types
        assert 'Xe~Xe' in ff.bond_types

        assert sympify('r') in ff.bond_types['Ar~Ar'].independent_variables
        assert ff.bond_types['Ar~Ar'].parameters['r_eq'] == u.unyt_quantity(10.0, u.nm)
        assert ff.bond_types['Ar~Ar'].parameters['k'] == u.unyt_quantity(10000, u.kJ / u.mol)
        assert ff.bond_types['Ar~Ar'].member_types == ('Ar', 'Ar')

        assert sympify('r') in ff.bond_types['Xe~Xe'].independent_variables
        assert ff.bond_types['Xe~Xe'].parameters['r_eq'] == u.unyt_quantity(10.0, u.nm)
        assert ff.bond_types['Xe~Xe'].parameters['k'] == u.unyt_quantity(20000, u.kJ / u.mol)
        assert ff.bond_types['Xe~Xe'].member_types == ('Xe', 'Xe')

    def test_ff_angletypes_from_xml(self, ff):
        assert len(ff.angle_types) == 2
        assert 'Ar~Ar~Ar' in ff.angle_types
        assert 'Xe~Xe~Xe' in ff.angle_types

        assert sympify('r') in ff.angle_types['Ar~Ar~Ar'].independent_variables
        assert ff.angle_types['Ar~Ar~Ar'].parameters['r_eq'] == u.unyt_quantity(10.0, u.nm)
        assert ff.angle_types['Ar~Ar~Ar'].parameters['z'] == u.unyt_quantity(100, u.kJ / u.mol)
        assert ff.angle_types['Ar~Ar~Ar'].member_types == ('Ar', 'Ar', 'Ar')

        assert sympify('r') in ff.angle_types['Xe~Xe~Xe'].independent_variables
        assert ff.angle_types['Xe~Xe~Xe'].parameters['r_eq'] == u.unyt_quantity(10.0, u.nm)
        assert ff.angle_types['Xe~Xe~Xe'].parameters['z'] == u.unyt_quantity(20, u.kJ / u.mol)
        assert ff.angle_types['Xe~Xe~Xe'].member_types == ('Xe', 'Xe', 'Xe')

    def test_ff_dihedraltypes_from_xml(self, ff):
        assert len(ff.dihedral_types) == 2
        assert 'Xe~Xe~Xe~Xe' in ff.dihedral_types
        assert 'Ar~Ar~Ar~Ar' in ff.dihedral_types

        assert sympify('r') in ff.dihedral_types['Ar~Ar~Ar~Ar'].independent_variables
        assert ff.dihedral_types['Ar~Ar~Ar~Ar'].parameters['r_eq'] == u.unyt_quantity(10.0, u.nm)
        assert ff.dihedral_types['Ar~Ar~Ar~Ar'].parameters['z'] == u.unyt_quantity(100, u.kJ / u.mol)
        assert ff.dihedral_types['Ar~Ar~Ar~Ar'].member_types == ('Ar', 'Ar', 'Ar', 'Ar')

        assert sympify('r') in ff.dihedral_types['Xe~Xe~Xe~Xe'].independent_variables
        assert ff.dihedral_types['Xe~Xe~Xe~Xe'].parameters['r_eq'] == u.unyt_quantity(10.0, u.nm)
        assert ff.dihedral_types['Xe~Xe~Xe~Xe'].parameters['z'] == u.unyt_quantity(20, u.kJ / u.mol)
        assert ff.dihedral_types['Xe~Xe~Xe~Xe'].member_types == ('Xe', 'Xe', 'Xe', 'Xe')

    def test_ff_impropertypes_from_xml(self, ff):
        assert len(ff.improper_types) == 1
        assert 'Xe~Xe~Xe~Xe' in ff.improper_types

        assert sympify('r') in ff.improper_types['Xe~Xe~Xe~Xe'].independent_variables
        assert ff.improper_types['Xe~Xe~Xe~Xe'].parameters['r_eq'] == u.unyt_quantity(10.0, u.nm)
        assert ff.improper_types['Xe~Xe~Xe~Xe'].parameters['z'] == u.unyt_quantity(20, u.kJ / u.mol)
        assert ff.improper_types['Xe~Xe~Xe~Xe'].member_types == ('Xe', 'Xe', 'Xe', 'Xe')


    def test_ff_charmm_xml(self):
        charm_ff = ForceField(get_path('trimmed_charmm.xml'))

        assert charm_ff.name == 'topologyCharmm'
        assert "*~CS~SS~*" in charm_ff.dihedral_types

        # Test list of parameters
        assert isinstance(charm_ff.dihedral_types["*~CE1~CE1~*"].parameters['k'], list)

        # This ensures that even though the parameters is a list, they can be hashed (by equality checks)
        assert charm_ff.dihedral_types["*~CE1~CE1~*"] == charm_ff.dihedral_types["*~CE1~CE1~*"]
        assert len(charm_ff.dihedral_types["*~CE1~CE1~*"].parameters['k']) == 2

        # Test Correct Parameter Values
        assert charm_ff.dihedral_types["*~CE1~CE1~*"].parameters['k'] == \
               [u.unyt_quantity(0.6276, u.kJ), u.unyt_quantity(35.564, u.kJ)]

    def test_non_unique_params(self):
        with pytest.raises(DocumentInvalid):
            ForceField(get_path('ff-example-nonunique-params.xml'))

    def test_missing_params(self):
        with pytest.raises(ForceFieldParseError):
            ForceField(get_path('ff-example-missing-parameter.xml'))

    def test_elementary_charge_to_coulomb(self, ff):
        elementary_charge = ff.atom_types['Li'].charge.to(u.elementary_charge)
        assert elementary_charge.units == u.Unit(u.elementary_charge)

    def test_atomclass_groups_charm_buck_ff(self):
        ff = ForceField(get_path('opls_charmm_buck.xml'))
        assert len(ff.atom_class_groups['CT']) == 2

    def test_ff_periodic_dihedrals_from_alphanumeric_symbols(self):
        ff = ForceField(get_path('opls_charmm_buck.xml'))
        assert 'A' in ff.atom_types['buck_O'].parameters
        with pytest.raises(TypeError):
            assert len(ff.dihedral_types['opls_140~*~*~opls_140'].parameters['c0'])
        assert len(ff.dihedral_types['NH2~CT1~C~O'].parameters['delta']) == 1

    def test_ff_from_etree(self):
        ff_etree = lxml.etree.parse(get_path('opls_charmm_buck.xml'))
        ff = ForceField(ff_etree)
        assert ff

    def test_ff_from_etree_iterable(self):
        ff_etrees = [
            lxml.etree.parse(get_path('opls_charmm_buck.xml')),
            lxml.etree.parse(get_path('trimmed_charmm.xml'))
        ]
        ff = ForceField(ff_etrees)
        assert ff

    def test_ff_mixed_type_error(self):
        with pytest.raises(TypeError):
            ff = ForceField([5, '20'])
            
    def test_named_potential_groups(self, named_groups_ff):
        assert named_groups_ff.potential_groups['BuckinghamPotential']
        assert named_groups_ff.angle_types['Xe~Xe~Xe'] in named_groups_ff.potential_groups['HarmonicAngle'].values()
        assert len(named_groups_ff.potential_groups['BuckinghamPotential']) == 3
        assert len(named_groups_ff.potential_groups['HarmonicBond']) == 2
        assert len(named_groups_ff.potential_groups['HarmonicAngle']) == 2
        assert len(named_groups_ff.potential_groups['PeriodicProper']) == 2
        assert len(named_groups_ff.potential_groups['RBProper']) == 1


    def test_potential_types_by_expression(self, named_groups_ff):
        atom_types_grouped_by_expression = named_groups_ff.group_atom_types_by_expression()
        bond_types_grouped_by_expression = named_groups_ff.group_bond_types_by_expression()
        angle_types_grouped_by_expression = named_groups_ff.group_angle_types_by_expression()
        dihedral_types_grouped_by_expression = named_groups_ff.group_dihedral_types_by_expression()
        improper_types_gropued_by_expression = named_groups_ff.group_improper_types_by_expression()

        assert len(atom_types_grouped_by_expression['A*exp(-B/r) - C/r**6']) == 3
        assert len(bond_types_grouped_by_expression['0.5*k*(r - r_eq)**2']) == 2
        assert len(angle_types_grouped_by_expression['0.5*z*(r - r_eq)**2']) == 2
        assert len(dihedral_types_grouped_by_expression['0.5*z*(r - r_eq)**2']) == 2
        assert len(improper_types_gropued_by_expression['0.5*z*(r - r_eq)**2']) == 1

