"""Calculate charges for a gmso Topology."""

from gmso.external.convert_openmm import to_openff_molecule
from gmso.utils.slicing import slice_topology_by_molecule


def calculate_charges(topology, method="am1bcc", slice_by="Molecule"):
    # Get all of the moleucles types in the topology
    charges_dict = dict()

    molecules = set()
    for site in topology.sites:
        molecules.add(site.molecule.name)

    # Figure out how to handle the tags
    for mol in molecules:
        mol_slice = slice_topology_by_molecule(topology=topology, molecule_tag=0)
        openff_mol = to_openff_molecule(mol_slice)
        openff_mol.assign_partial_charges(partial_charge_method=method)
        for atom in openff_mol.atoms:
            print(atom.partial_charge)
