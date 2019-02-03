import mbuild as mb
from mbuild.examples import Ethane

from topology.external.convert_mbuild import from_mbuild, to_mbuild

def test_from_mbuild():
    ethane = Ethane()
    top = from_mbuild(ethane)

    assert top.n_sites == 8
    assert top.n_connections == 7

def test_full_conversion():
    ethane = Ethane()
    top = from_mbuild(ethane)

    new = to_mbuild(top)

    assert new.n_particles == 8
    assert new.n_bonds == 7
