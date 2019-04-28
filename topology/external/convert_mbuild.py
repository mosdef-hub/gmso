import unyt as u

from topology.core.topology import Topology
from topology.core.subtopology import SubTopology
from topology.core.site import Site
from topology.core.bond import Bond
from topology.utils.io import has_mbuild


if has_mbuild:
    import mbuild as mb

def from_mbuild(compound):
    msg = ("Provided argument that is not an mbuild Compound")
    assert isinstance(compound, mb.Compound), msg

    top = Topology()

    # Keep the name if it is not the default mbuild Compound name
    if compound.name != mb.Compound().name:
        top.name = compound.name

    site_map = dict()
    for child in compound.children:
        if child.children is None:
            pos = child.xyz[0] * u.nanometer
            site = Site(name=child.name, position=pos)
            site_map[child] = site
        else:
            subtop = SubTopology(name=child.name)
            top.add_subtopology(subtop)
            for childchild in child.children:
                pos = childchild.xyz[0] * u.nanometer
                site = Site(name=childchild.name, position=pos)
                site_map[childchild] = site
                subtop.add_site(site)
    top.update_top()

    for b1, b2 in compound.bonds():
        new_bond = Bond(connection_members=[site_map[b1], site_map[b2]],
                connection_type=None)
        top.add_connection(new_bond, update_types=False)
    top.update_top()

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
    for site in topology.sites:
        particle = mb.Compound(name=site.name, pos=site.position[0])
        particle_map[site] = particle
        compound.add(particle)

    for connect in topology.connections:
        if isinstance(connect, Bond):
            compound.add_bond((
                particle_map[connect.connection_members[0]],
                particle_map[connect.connection_members[1]]))

    return compound
