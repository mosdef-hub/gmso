

from topology.forcefield import ForceField
from topology.tests.utils import get_path
from topology.tests.base_test import BaseTest

class TestForceField(BaseTest):

    def test_ff_atomtypes(self):
        ForceField.from_xml(get_path('new_file.xml'))

