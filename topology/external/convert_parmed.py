import parmed as pmd
#import unyt as u

from topology.core.topology import Topology
from topology.core.site import Site
#from topology.core.box import Box

def from_parmed(structure):
    msg = ("Provided argument that is not a Parmed Structure")
    assert isinstance(structure, pmd.Structure), msg

    top = Topology(name=structure.title)
    map = dict()
    for atom in structure.atoms:
        site = Site(name=atom.name, position=[atom.xx, atom.xy, atom.xz])
        #site = Site(name=atom.name, position=[atom.xx, atom.xy, atom.xz]*u.nanometer)
        map[atom] = site
        top.add_site(site)

    #if structure.box:
        # top.box = Box(structure.box[0:3]*u.nanometer, angles=structure.box[4:7])

    for bond in structure.bonds:
        if map[bond.atom2] not in map[bond.atom1].connections:
            map[bond.atom1].add_connection(map[bond.atom2])
        if map[bond.atom1] not in map[bond.atom2].connections:
            map[bond.atom2].add_connection(map[bond.atom1])

    top.update_connection_list()

    return top

