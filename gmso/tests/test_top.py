import pytest

import gmso
import parmed as pmd
from gmso.formats.top import write_top
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn
from gmso.tests.utils import get_path
from gmso.exceptions import EngineIncompatibilityError
from gmso.external import from_mbuild


class TestTop(BaseTest):
    def test_write_top(self, typed_ar_system):
        top = typed_ar_system
        write_top(top, 'ar.top')


    @pytest.mark.parametrize('top', ['typed_ar_system',
        'typed_water_system'])
    def test_pmd_loop(self, top, request):
        write_top(request.getfixturevalue(top), 'system.top')
        pmd.load_file('system.top')


    def test_modified_potentials(self, ar_system):
        top = from_mbuild(ar_system())

        ff = gmso.ForceField(get_fn('ar.xml'))

        for site in top.sites:
            site.atom_type = ff.atom_types['Ar']

        top.update_topology()

        top.atom_types[0].set_expression('sigma + epsilon')

        with pytest.raises(EngineIncompatibilityError):
            write_top(top, 'out.top')

        alternate_lj = '4*epsilon*sigma**12/r**12 - 4*epsilon*sigma**6/r**6'
        top.atom_types[0].set_expression(alternate_lj)

        write_top(top, 'ar.top')


    def test_water_top(self, water_system):
        top = water_system

        ff = gmso.ForceField(get_path('tip3p.xml'))

        for site in top.sites:
            site.atom_type = ff.atom_types[site.name]

        top.update_sites()
        top.update_atom_types()

        for bond in top.bonds:
            bond.bond_type = bond.connection_type = ff.bond_types['opls_111~opls_112']

        top.update_bonds()
        top.update_bond_types()

        for subtop in top.subtops:
            angle = gmso.core.angle.Angle(
                connection_members=[site for site in subtop.sites],
                name="opls_112~opls_111~opls_112",
                connection_type=ff.angle_types["opls_112~opls_111~opls_112"]
            )
            top.add_connection(angle)

        top.update_angles()
        top.update_angle_types()

        write_top(top, 'water.top')
