import pytest

import topology as topo
from topology.formats.top import write_top
from topology.tests.base_test import BaseTest
from topology.utils.io import get_fn
from topology.exceptions import EngineIncompatibilityError


class TestTop(BaseTest):
    def test_write_top(self, ar_system):
        top = ar_system

        ff = topo.ForceField(get_fn('ar.xml'))

        for site in top.sites:
            site.atom_type = ff.atom_types['Ar']

        top.update_topology()

        write_top(top, 'ar.top')


    def test_modified_potentials(self, ar_system):
        top = ar_system

        ff = topo.ForceField(get_fn('ar.xml'))

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

        ff = topo.ForceField(get_fn('topology-tip3p.xml'))

        for site in top.sites:
            site.atom_type = ff.atom_types[site.name]

        top.dupate_topology()

        write_top(top, 'water.top')
