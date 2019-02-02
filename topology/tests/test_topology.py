from topology.core.topology import Topology
from topology.core.site import Site

def test_new_topology():
    top = Topology(name='mytop')
    assert top.name == 'mytop'

def test_add_site():
    top = Topology()
    site = Site(name='site')

    assert top.n_sites == 0
    top.add_site(site)
    assert top.n_sites == 1
