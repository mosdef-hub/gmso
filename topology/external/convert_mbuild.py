import mbuild as mb

from topology.core.topology import Topology
from topology.core.site import Site

def from_mbuild(compound):
    msg = ("Provided argument that is not an mbuild Compound")
    assert isinstance(compound, mb.Compound), msg

    top = Topology(name=compound.name)
    map = dict()
    for child in compound.particles():
        site = Site(name=child.name, position=child.xyz)
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

def to_mbuild(topology):
    msg = ("Provided argument that is not a topology")
    assert isinstance(topology, Topology), msg

    compound = mb.Compound()
    if topology.name is None:
        compound.name = 'Compound'
    else:
        compound.name = topology.name

    map = dict()
    for site in topology.site_list:
        particle = mb.Compound(name=site.name, pos=site.position[0])
        map[site] = particle
        compound.add(particle)

    for connect in topology.connection_list:
        compound.add_bond((map[connect.site1], map[connect.site2]))

    return compound
