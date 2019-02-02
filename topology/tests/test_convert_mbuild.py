import mbuild as mb
from mbuild.examples import Ethane

from topology.external.convert_mbuild import from_mbuild

def test_from_mbuild():
    ethane = Ethane()
    top = from_mbuild(ethane)
    import pdb; pdb.set_trace()

    assert top.n_sites == 8
