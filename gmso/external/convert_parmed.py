import numpy as np
import unyt as u
import sympy as sym
from sympy.parsing.sympy_parser import parse_expr
import warnings

import gmso
from gmso.utils.io import import_, has_parmed
from topology.core.element import element_by_name, element_by_symbol, element_by_atom_type

if has_parmed:
    pmd = import_('parmed')

def from_parmed(structure):
    """Convert a parmed.Structure to a gmso.Topology

    Parameters
    ----------
    structure : parmed.Structure
        parmed.Structure instance that need to be converted

    Returns
    -------
    top : gmso.Topology
    """
    msg = ("Provided argument is not a Parmed Structure")
    assert isinstance(structure, pmd.Structure), msg

    top = gmso.Topology(name=structure.title)
    site_map = dict()
    # simple check if our pmd.Structure is fully parametrized
    is_parametrized = True if isinstance(structure.atoms[0].atom_type,
            pmd.AtomType) else False

    # Consolidate parmed atomtypes and relate topology atomtypes
    if is_parametrized:
        unique_atom_types = list(set([a.atom_type for a in structure.atoms]))
        pmd_top_atomtypes = {}
        for atom_type in unique_atom_types:
            top_atomtype = topo.AtomType(
                    name=atom_type.name,
                    charge=atom_type.charge * u.elementary_charge,
                    parameters={
                        'sigma': (atom_type.sigma * u.angstrom).in_units(u.nm),
                        'epsilon': atom_type.epsilon * u.Unit('kcal / mol')
                        })
            pmd_top_atomtypes[atom_type] = top_atomtype

        # Consolidate parmed bondtypes and relate to topology bondtypes
        pmd_top_bondtypes = {}
        for btype in structure.bond_types:
            bond_params = {
                    'k': (2 * btype.k * u.Unit('kcal / (nm**2 * mol)')),
                    'r_eq': (btype.req * u.angstrom).in_units(u.nm)
                    }

            top_bondtype = topo.BondType(parameters=bond_params)
            pmd_top_bondtypes[btype] = top_bondtype

        # Consolidate parmed angletypes and relate to topology angletypes
        pmd_top_angletypes = {}
        for angletype in structure.angle_types:
            angle_params = {
                    'k': (2 * angletype.k * u.Unit('kcal / (rad**2 * mol)')),
                    'theta_eq': (angletype.theteq * u.degree)
                }
                 
            top_angletype = topo.AngleType(parameters=angle_params)
            pmd_top_angletypes[angletype] = top_angletype

        # Consolidate parmed dihedraltypes and relate to topology dihedraltypes
        pmd_top_dihedraltypes = {}
        for dihedraltype in structure.dihedral_types:
            dihedral_params = {
                    'k': (dihedraltype.phi_k * u.Unit('kcal / mol')),
                    'phi_eq': (dihedraltype.phase * u.degree),
                    'n': dihedraltype.per * u.dimensionless
                }
                 
            top_dihedraltype = topo.DihedralType(parameters=dihedral_params)
            pmd_top_dihedraltypes[dihedraltype] = top_dihedraltype
        for dihedraltype in structure.rb_torsion_types:
            dihedral_params = {
                    'c0': (dihedraltype.c0 * u.Unit('kcal/mol')),
                    'c1': (dihedraltype.c1 * u.Unit('kcal/mol')),
                    'c2': (dihedraltype.c2 * u.Unit('kcal/mol')),
                    'c3': (dihedraltype.c3 * u.Unit('kcal/mol')),
                    'c4': (dihedraltype.c4 * u.Unit('kcal/mol')),
                    'c5': (dihedraltype.c5 * u.Unit('kcal/mol')),
                }
                 
            top_dihedraltype = topo.DihedralType(parameters=dihedral_params,
                    expression='c0 * cos(phi)**0 + c1 * cos(phi)**1 + ' +
                    'c2 * cos(phi)**2 + c3 * cos(phi)**3 + c4 * cos(phi)**4 + ' +
                    'c5 * cos(phi)**5',
                    independent_variables='phi'
                    )
            pmd_top_dihedraltypes[dihedraltype] = top_dihedraltype

    for atom in structure.atoms:
        if is_parametrized:
            site = topo.Site(
                name=atom.name,
                charge=atom.charge * u.elementary_charge,
                position=([atom.xx, atom.xy, atom.xz] * u.angstrom).in_units(
                    u.nm),
                atom_type=pmd_top_atomtypes[atom.atom_type])
        else:
            site = gmso.Site(
                name=atom.name,
                charge=atom.charge * u.elementary_charge,
                position=([atom.xx, atom.xy, atom.xz] * u.angstrom).in_units(
                    u.nm),
                atom_type=None)
        site_map[atom] = site
        top.add_site(site)
    top.update_topology()

    if np.all(structure.box):
        # This is if we choose for topology to have abox
        top.box = gmso.Box(
            (structure.box[0:3] * u.angstrom).in_units(u.nm),
            angles=u.degree * structure.box[3:6])

    for bond in structure.bonds:
        # Generate bond parameters for BondType that gets passed
        # to Bond
        if is_parametrized:
            top_connection = gmso.Bond(connection_members=[site_map[bond.atom1],

                site_map[bond.atom2]],
                connection_type=pmd_top_bondtypes[bond.type])

        # No bond parameters, make Connection with no connection_type
        else:
            top_connection = gmso.Bond(connection_members=[site_map[bond.atom1],
                site_map[bond.atom2]],
                connection_type=None)

        top.add_connection(top_connection, update_types=False)
    top.update_topology()

    for angle in structure.angles:
        # Generate angle parameters for AngleType that gets passed
        # to Angle
        if is_parametrized:
            top_connection = gmso.Angle(connection_members=[site_map[angle.atom1],
                site_map[angle.atom2], site_map[angle.atom3]],
                connection_type=pmd_top_angletypes[angle.type])

        # No bond parameters, make Connection with no connection_type
        else:
            top_connection = gmso.Angle(connection_members=[site_map[angle.atom1],
                site_map[angle.atom2], site_map[angle.atom3]],
                connection_type=None)

        top.add_connection(top_connection, update_types=False)
    top.update_topology()

    for dihedral in structure.dihedrals:
        # Generate dihedral parameters for DihedralType that gets passed
        # to Dihedral
        # These all follow periodic torsions functions
        # (even if they are improper dihedrals)
        # Which are the default expression in top.DihedralType
        # These periodic torsion dihedrals get stored in top.dihedrals
        if dihedral.improper:
            warnings.warn("ParmEd improper dihedral {} ".format(dihedral) +
                    "following periodic torsion " +
                    "expression detected, currently accounted for as " +
                    "topology.Dihedral with a periodic torsion expression")
        if is_parametrized:
            top_connection = topo.Dihedral(connection_members=
                    [site_map[dihedral.atom1], site_map[dihedral.atom2], 
                        site_map[dihedral.atom3], site_map[dihedral.atom4]],
                    connection_type=pmd_top_dihedraltypes[dihedral.type])

        # No bond parameters, make Connection with no connection_type
        else:
            top_connection = topo.Dihedral(connection_members=
                    [site_map[dihedral.atom1], site_map[dihedral.atom2], 
                        site_map[dihedral.atom3], site_map[dihedral.atom4]],
                connection_type=None)

        top.add_connection(top_connection, update_types=False)

    for rb_torsion in structure.rb_torsions:
        # Generate dihedral parameters for DihedralType that gets passed
        # to Dihedral
        # These all follow RB torsion functions
        # These RB torsion dihedrals get stored in top.dihedrals
        if rb_torsion.improper:
            warnings.warn("ParmEd improper dihedral {} ".format(rb_torsion) +
                    "following RB torsion " +
                    "expression detected, currently accounted for as " +
                    "topology.Dihedral with a RB torsion expression")
        if is_parametrized:
            top_connection = topo.Dihedral(connection_members=
                    [site_map[rb_torsion.atom1], site_map[rb_torsion.atom2], 
                        site_map[rb_torsion.atom3], site_map[rb_torsion.atom4]],
                    connection_type=pmd_top_dihedraltypes[rb_torsion.type])

        # No bond parameters, make Connection with no connection_type
        else:
            top_connection = topo.Dihedral(connection_members=
                    [site_map[rb_torsion.atom1], site_map[rb_torsion.atom2], 
                        site_map[rb_torsion.atom3], site_map[rb_torsion.atom4]],
                connection_type=None)

        top.add_connection(top_connection, update_types=False)

    top.update_topology()

    return top


def to_parmed(top):
    """Convert a topology.Topology to a parmed.Structure

    At this point we only assume a three level structure for topology
    Topology - Subtopology - Sites, which transform to three level of
    Parmed Structure - Residue - Atoms.
    If we decide to support multiple level Subtopology in the future,
    this method will need some re-work. Tentative plan is to have the
    Parmed Residue to be equivalent to the Subtopology right above Site.
    Parameters
    ----------
    top : topology.Topology
        topology.Topology instance that need to be converted

    Returns
    -------
    structure : parmed.Structure
    """
    # Sanity check
    msg = "Provided argument is not a topology.Topology."
    assert isinstance(top, topo.Topology)

    # Set up Parmed structure and define general properties
    structure = pmd.Structure()
    structure.title = top.name
    structure.box = np.concatenate((top.box.lengths.to('angstrom').value,
                                    top.box.angles.to('degree').value))

    # Set up unparametrized system
    # Build residue-atom map
    residue_map = dict()
    default_residue = pmd.Residue('RES')
    for subtop in top.subtops:
        for site in subtop.sites:
            residue_map[site] = subtop

    # Build up atom
    atom_mapping = dict()
    for site in top.sites:
        if site in residue_map:
            residue = residue_map[site]
        else:
            residue = default_residue
        # Check element
        if site.element:
            atomic_number = site.element.atomic_number
            charge = site.element.charge
        else:
            atomic_number = 0
            charge = 0

        pmd_atom = pmd.Atom(atomic_number=atomic_number, name=site.name,
                            mass=site.mass, charge=site.charge)
        pmd_atom.xx, pmd_atom.xy, pmd_atom.xz = site.position.to('angstrom').value

        # Add atom to structure
        structure.add_atom(pmd_atom, resname=residue.name, resnum=residue.idx)
        atom_mapping[atom] = pmd_atom

    # "Claim" all of the item it contains and subsequently index all of its item
    structure.residues.claim()

    # Create and add bonds to Parmed structure
    for bond in top.bonds:
        atom1, atom2 = bond.connection_members
        bond = pmd.Bond(atom_mapping[atom1], atom_mapping[atom2])
        structure.bonds.append(bond)

    # Set up structure for Connection Type conversion
    # Will need to work on these helper function later
    if top.atom_types:
        _convert_atom_types(top, structure)
    if top.bond_types:
        _convert_bond_types(top, structure)
    if top.angle_types:
        _convert_angle_types(top, structure)
    if top.dihedral_types:
        _convert_dihedral_types(top, structure)

    return structure


def _convert_atom_types(top, structure):
    """Helper function to convert Topology AtomType to Structure AtomType

    This function will first check the AtomType expression of Topology and make sure it match with the one default in Parmed.
    After that, it would start atomtyping and parametrizing this part of the structure.

    Parameters
    ----------
    top : topology.Topology
        The topology that need to be converted
    structure: parmed.Structure
        The destination parmed Structure
    """
    pass


def _convert_bond_types(top, structure):
    """Helper function to convert Topology BondType to Structure BondType

    This function will first check the BondType expression of Topology and make sure it match with the one default in Parmed.
    After that, it would start atomtyping and parametrizing this part of the structure.

    Parameters
    ----------
    top : topology.Topology
        The topology that need to be converted
    structure: parmed.Structure
        The destination parmed Structure
    """
    pass


def _convert_atom_types(top, structure):
    """Helper function to convert Topology AngleType to Structure AngleType

    This function will first check the AngleType expression of Topology and make sure it match with the one default in Parmed.
    After that, it would start atomtyping and parametrizing the structure.

    Parameters
    ----------
    top : topology.Topology
        The topology that need to be converted
    structure: parmed.Structure
        The destination parmed Structure
    """
    pass


def _convert_atom_types(top, structure):
    """Helper function to convert Topology DihedralType to Structure DihedralType

    This function will first check the DihedralType expression of Topology and make sure it match with the one default in Parmed.
    After that, it would start atomtyping and parametrizing the structure.

    Parameters
    ----------
    top : topology.Topology
        The topology that need to be converted
    structure: parmed.Structure
        The destination parmed Structure
    """
    pass
