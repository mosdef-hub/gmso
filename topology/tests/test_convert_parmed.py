import parmed as pmd

from topology.external.convert_parmed import from_parmed
from topology.tests.base_test import BaseTest
from topology.utils.io import get_fn


class TestConvertParmEd(BaseTest):

    def test_from_parmed_basic(self):
        struc = pmd.load_file(get_fn('ethane.mol2')).to_structure()
        top = from_parmed(struc)
        for site in top.site_list:
            assert site.atom_type is None 
        for connection in top.connection_list:
            assert connection.connection_type is None 
        assert top.n_sites == 8
        assert top.n_connections == 7
    
    def test_from_parmed_parametrized_structure(self):
        struc = pmd.load_file(get_fn('ethane.top'), xyz=get_fn('ethane.gro'))
        top = from_parmed(struc)
        assert top.n_sites == 8
        assert top.n_connections == 7
    
        for site in top.site_list:
            assert site.atom_type is not None
            assert site.charge is not None
    
        for connection in top.connection_list:
            assert connection.connection_type is not None 
