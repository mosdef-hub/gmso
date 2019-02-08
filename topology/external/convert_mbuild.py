import mbuild as mb

from topology.core.topology import Topology
from topology.core.site import Site

def from_mbuild(compound):
    msg = ("Provided argument that is not an mbuild Compound")
    assert isinstance(compound, mb.Compound), msg

    top = Topology(name=compound.name)
    site_map = dict()
    for child in compound.particles():
        site = Site(name=child.name, position=child.xyz)
        site_map[child] = site
        top.add_site(site)

    for b1, b2 in compound.bonds():
        if site_map[b2] not in site_map[b1].connections:
            site_map[b1].add_connection(site_map[b2])
        if site_map[b1] not in site_map[b2].connections:
            site_map[b2].add_connection(site_map[b1])

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

    particle_map = dict()
    for site in topology.site_list:
        particle = mb.Compound(name=site.name, pos=site.position[0])
        particle_map[site] = particle
        compound.add(particle)

    for connect in topology.connection_list:
        compound.add_bond((particle_map[connect.site1], particle_map[connect.site2]))

    return compound
