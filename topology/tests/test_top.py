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

        top.atom_types[0].expression = 'sigma + epsilon'

        with pytest.raises(EngineIncompatibilityError):
            write_top(top, 'out.top')
