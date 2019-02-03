import mbuild as mb
from simtk.unit import nanometer

from topology.core.topology import Topology
from topology.core.site import Site

def from_mbuild(compound):
    msg = ("Provided argument that is not an mbuild Compound")
    assert isinstance(compound, mb.Compound), msg

    top = Topology(name=compound.name)
    map = dict()
    for child in compound.particles():
        pos = [val * nanometer for val in child.xyz]
        site = Site(name=child.name, position=pos)
        map[child] = site
        top.add_site(site)

    for b1, b2 in compound.bonds():
        #top.add_bond(map[b1], map[b1])
        if map[b2] not in map[b1].connections:
            map[b1].add_connection(map[b2])
        if map[b1] not in map[b2].connections:
            map[b2].add_connection(map[b1])

    top.update_connection_list()

    return top
