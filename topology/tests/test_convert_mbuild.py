import mbuild as mb
from mbuild.examples import Ethane

from topology.external.convert_mbuild import from_mbuild

def test_from_mbuild():
    ethane = Ethane()
    top = from_mbuild(ethane)

    assert top.n_sites == 8
    assert top.n_connections == 7
