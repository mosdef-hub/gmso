

from topology.forcefield import ForceField
from topology.tests.utils import get_path
from topology.tests.base_test import BaseTest


class TestForceField(BaseTest):

    def test_ff_atomtypes(self):
        ff = ForceField.from_xml(get_path('matt-at.xml'))
        assert len(ff.atom_types) == 2
