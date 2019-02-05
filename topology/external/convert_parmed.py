import parmed as pmd

from topology.core.topology import Topology
from topology.core.site import Site

def from_parmed(structure):
    msg = ("Provided argument that is not a Parmed Structure")
    assert isinstance(structure, pmd.Structure), msg

    top = Topology(name=structure.name)
    map = dict()
    for atom in structure.atoms:
        site = Site(name=atom.name, position=[atom.xx, atom.xy, atom.xz])
        map[atom] = site
        top.add_site(site)

    for bond in structure.bonds:
        if map[bond.atom2] not in map[bond.atom1].connections:
            map[bond.atom1].add_connection(map[bond.atom2])
        if map[bond.atom1] not in map[bond.atom2].connections:
            map[bond.atom2].add_connection(map[bond.atom1])

    top.update_connection_list()

    return top

