import pytest

import topology as topo
from topology.formats.mcf import write_mcf
from topology.tests.base_test import BaseTest
from topology.utils.io import get_fn
from topology.exceptions import EngineIncompatibilityError
from topology.external.convert_mbuild import from_mbuild

class TestTop(BaseTest):
    def test_write_lj_mcf(self, single_ar):
        top = single_ar

        ff = topo.ForceField(get_fn('ar.xml'))

        for site in top.sites:
            site.atom_type = ff.atom_types['Ar']

        top.update_topology()

        write_mcf(top, 'ar.mcf')

    def test_write_mie_mcf(self, single_xe):
        top = single_xe

        ff = topo.ForceField(get_fn('xe_mie.xml'))

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
        assert mcf_data[atom_section_start+2][3] == ''
        assert mcf_data[atom_section_start+2][4] == '0.0'
        assert mcf_data[atom_section_start+2][5] == 'Mie'
        assert np.isclose(float(mcf_data[atom_section_start+2][6]),)
        assert np.isclose(float(mcf_data[atom_section_start+2][7]),)
        assert np.isclose(float(mcf_data[atom_section_start+2][8]),)
        assert np.isclose(float(mcf_data[atom_section_start+2][9]),)

 
    def test_modified_potentials(self, single_ar):
        top = single_ar

        ff = topo.ForceField(get_fn('ar.xml'))

        for site in top.sites:
            site.atom_type = ff.atom_types['Ar']

        top.update_topology()

        top.atom_types[0].set_expression('sigma + epsilon')

        with pytest.raises(EngineIncompatibilityError):
            write_mcf(top, 'out.mcf')

        alternate_lj = '4*epsilon*sigma**12/r**12 - 4*epsilon*sigma**6/r**6'
        top.atom_types[0].set_expression(alternate_lj)

        write_mcf(top, 'ar.mcf')
