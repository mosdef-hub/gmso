import numpy as np
import parmed as pmd
import unyt as u

from topology.core.topology import Topology
from topology.core.site import Site
from topology.core.atom_type import AtomType
from topology.core.bond import Bond
from topology.core.bond_type import BondType
from topology.core.box import Box
from topology.tests.base_test import BaseTest


def from_parmed(structure):
    msg = ("Provided argument that is not a Parmed Structure")
    assert isinstance(structure, pmd.Structure), msg

    top = Topology(name=structure.title)
    site_map = dict()
    for atom in structure.atoms:
        if isinstance(atom.atom_type, pmd.AtomType):
            atom_type = AtomType(
                name=atom.atom_type.name,
                charge=atom.atom_type.charge * u.elementary_charge,
                parameters={
                    'sigma': (atom.sigma * u.angstrom).in_units(u.nm),
                    'epsilon': atom.epsilon * u.Unit('kcal / mol')
                })
            site = Site(
                name=atom.name,
                charge=atom.charge * u.elementary_charge,
                position=([atom.xx, atom.xy, atom.xz] * u.angstrom).in_units(
                    u.nm),
                atom_type=atom_type)
        else:
            site = Site(
                name=atom.name,
                charge=atom.charge * u.elementary_charge,
                position=([atom.xx, atom.xy, atom.xz] * u.angstrom).in_units(
                    u.nm),
                atom_type=None)
        site_map[atom] = site
        top.add_site(site)

    if np.all(structure.box):
        # This is if we choose for topology to have abox
        top.box = Box(
            (structure.box[0:3] * u.angstrom).in_units(u.nm),
            angles=u.degree * structure.box[3:6])

    for bond in structure.bonds:
        # Generate bond parameters for BondType that gets passed
        # to Bond
        if isinstance(bond.type, pmd.BondType):
            bond_params = {
                'k': (2 * bond.type.k * u.Unit('kcal / (angstrom**2 * mol)')),
                'r_eq': (bond.type.req * u.angstrom).in_units(u.nm)
            }
            new_connection_type = BondType(parameters=bond_params)
            top_connection = Bond(bond_partners=[site_map[bond.atom1], 
                site_map[bond.atom2]],
                connection_type=new_connection_type)

        # No bond parameters, make Connection with no connection_type
        else:
            top_connection = Bond(bond_parnters=[site_map[bond.atom1],
                site_map[bond.atom2]],
                connection_type=None)

        top.add_connection(top_connection)

        #if site_map[bond.atom2] not in site_map[bond.atom1].connections:
        #    site_map[bond.atom1].add_connection(site_map[bond.atom2])
        #if site_map[bond.atom1] not in site_map[bond.atom2].connections:
        #    site_map[bond.atom2].add_connection(site_map[bond.atom1])

    # TODO: Angles
    # TODO: Dihedrals

    top.update_connection_list()

    return top
