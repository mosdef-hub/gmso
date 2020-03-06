import pytest
import mbuild
import numpy as np
import unyt as u

import gmso 
from gmso.formats.mcf import write_mcf
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path
from gmso.utils.io import get_fn
from gmso.exceptions import EngineIncompatibilityError
from gmso.external.convert_mbuild import from_mbuild

class TestMCF(BaseTest):
    def test_write_lj_mcf(self):
        top = from_mbuild(mbuild.Compound(name='Ar'))

        ff = gmso.ForceField(get_fn('ar.xml'))

        for site in top.sites:
            site.atom_type = ff.atom_types['Ar']

        top.update_topology()

        write_mcf(top, 'ar.mcf')

        mcf_data = []
        with open('ar.mcf') as f:
            for line in f:
                mcf_data.append(line.strip().split())

        for idx,line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == 'Atom_Info':
                    atom_section_start = idx
                elif line[1] == 'Bond_Info':
                    bond_section_start = idx
                elif line[1] == 'Angle_Info':
                    angle_section_start = idx
                elif line[1] == 'Dihedral_Info':
                    dihedral_section_start = idx
                elif line[1] == 'Improper_Info':
                    improper_section_start = idx
                elif line[1] == 'Fragment_Info':
                    fragment_section_start = idx
                elif line[1] == 'Fragment_Connectivity':
                    fragment_conn_start = idx

        assert mcf_data[atom_section_start+1][0] == '1'
        assert mcf_data[atom_section_start+2][1] == 'Ar'
        assert mcf_data[atom_section_start+2][2] == 'Ar'
        assert mcf_data[atom_section_start+2][5] == 'LJ'
        assert np.isclose(float(mcf_data[atom_section_start+2][3]),
                          ff.atom_types['Ar'].mass.in_units(u.amu).value)
        assert np.isclose(float(mcf_data[atom_section_start+2][4]),
                          ff.atom_types['Ar'].charge.in_units(u.elementary_charge).value)
        assert np.isclose(float(mcf_data[atom_section_start+2][6]),
                          (ff.atom_types['Ar'].parameters['epsilon'] / u.kb).in_units(u.K).value)
        assert np.isclose(float(mcf_data[atom_section_start+2][7]),
                          ff.atom_types['Ar'].parameters['sigma'].in_units(u.Angstrom).value)


    def test_write_mie_mcf(self):
        top = from_mbuild(mbuild.Compound(name='Xe'))

        ff = gmso.ForceField(get_path('noble_mie.xml'))

        for site in top.sites:
            site.atom_type = ff.atom_types['Xe']

        top.update_topology()

        write_mcf(top, 'xe.mcf')

        mcf_data = []
        with open('xe.mcf') as f:
            for line in f:
                mcf_data.append(line.strip().split())

        for idx,line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == 'Atom_Info':
                    atom_section_start = idx
                elif line[1] == 'Bond_Info':
                    bond_section_start = idx
                elif line[1] == 'Angle_Info':
                    angle_section_start = idx
                elif line[1] == 'Dihedral_Info':
                    dihedral_section_start = idx
                elif line[1] == 'Improper_Info':
                    improper_section_start = idx
                elif line[1] == 'Fragment_Info':
                    fragment_section_start = idx
                elif line[1] == 'Fragment_Connectivity':
                    fragment_conn_start = idx

        # Check a some atom info
        assert mcf_data[atom_section_start+1][0] == '1'
        assert mcf_data[atom_section_start+2][1] == 'Xe'
        assert mcf_data[atom_section_start+2][2] == 'Xe'
        assert mcf_data[atom_section_start+2][5] == 'Mie'
        assert np.isclose(float(mcf_data[atom_section_start+2][3]),
                          ff.atom_types['Xe'].mass.in_units(u.amu).value)
        assert np.isclose(float(mcf_data[atom_section_start+2][4]),
                          ff.atom_types['Xe'].charge.in_units(u.elementary_charge).value)
        assert np.isclose(float(mcf_data[atom_section_start+2][6]),
                          (ff.atom_types['Xe'].parameters['epsilon'] / u.kb).in_units(u.K).value)
        assert np.isclose(float(mcf_data[atom_section_start+2][7]),
                          ff.atom_types['Xe'].parameters['sigma'].in_units(u.Angstrom).value)
        assert np.isclose(float(mcf_data[atom_section_start+2][8]),
                          ff.atom_types['Xe'].parameters['n'])
        assert np.isclose(float(mcf_data[atom_section_start+2][9]),
                          ff.atom_types['Xe'].parameters['m'])

 
    def test_modified_potentials(self):
        top = from_mbuild(mbuild.Compound(name='Ar'))

        ff = gmso.ForceField(get_fn('ar.xml'))

        for site in top.sites:
            site.atom_type = ff.atom_types['Ar']

        top.update_topology()

        top.atom_types[0].set_expression('sigma + epsilon')

        with pytest.raises(EngineIncompatibilityError):
            write_mcf(top, 'out.mcf')

        alternate_lj = '4*epsilon*sigma**12/r**12 - 4*epsilon*sigma**6/r**6'
        top.atom_types[0].set_expression(alternate_lj)

        write_mcf(top, 'ar.mcf')
