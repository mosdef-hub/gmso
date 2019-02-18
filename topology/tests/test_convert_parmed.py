import parmed as pmd

from topology.external.convert_parmed import from_parmed
from topology.tests.base_test import BaseTest
from topology.utils.io import get_fn


class TestConvertParmEd(BaseTest):

    def test_from_parmed(self):
        struc = pmd.load_file(get_fn('ethane.mol2')).to_structure()
        top = from_parmed(struc)
    
        assert top.n_sites == 8
        assert top.n_connections == 7
