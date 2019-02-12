import parmed as pmd

from topology.external.convert_parmed import from_parmed
from topology.utils.io import get_fn


def test_from_parmed_basic():
    struc = pmd.load_file(get_fn('ethane.mol2')).to_structure()
    top = from_parmed(struc)

    assert top.n_sites == 8
    assert top.n_connections == 7

def test_from_parmed_parametrized_structure():
    struc = pmd.load_file(get_fn('ethane.top'), xyz=get_fn('ethane.gro'))
    top = from_parmed(struc)

    # Some tests for charges and atom_types for each site
