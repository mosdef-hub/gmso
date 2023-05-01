"""Module support for converting to/from ParmEd objects."""
import warnings
from operator import attrgetter

import numpy as np
import unyt as u
from sympy.parsing.sympy_parser import parse_expr

import gmso
from gmso.core.element import element_by_atom_type, element_by_atomic_number
from gmso.exceptions import GMSOError
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.utils.io import has_parmed, import_

if has_parmed:
    pmd = import_("parmed")

lib = PotentialTemplateLibrary()
# function to check reversibility of dihedral type
rev_dih_order = lambda x: x.atom2.type > x.atom3.type or (
    x.atom2.type == x.atom3.type and x.atom1.type > x.atom4.type
)
# function to get the names for a given atomtype
atomtype_namegetter = lambda x: (
    attrgetter("atom1.type", "atom2.type", "atom3.type", "atom4.type")(x)
)


def from_parmed(structure, refer_type=True):
    """Convert a parmed.Structure to a gmso.Topology.

    Convert a parametrized or un-parametrized parmed.Structure object to a topology.Topology.
    Specifically, this method maps Structure to Topology and Atom to Site.
    This method can only convert AtomType, BondType AngleType, DihedralType, and
    ImproperType.

    Parameters
    ----------
    structure : parmed.Structure
        parmed.Structure instance that need to be converted.
    refer_type : bool, optional, default=True
        Whether or not to transfer AtomType, BondType, AngleType,
        DihedralType, and ImproperType information

    Returns
    -------
    top : gmso.Topology
    """
    msg = "Provided argument is not a Parmed Structure"
    assert isinstance(structure, pmd.Structure), msg

    top = gmso.Topology(name=structure.title)
    site_map = dict()

    if np.all(structure.box):
        # This is if we choose for topology to have abox
        top.box = gmso.Box(
            (structure.box[0:3] * u.angstrom).in_units(u.nm),
            angles=u.degree * structure.box[3:6],
        )

    # Consolidate parmed atomtypes and relate topology atomtypes
    if refer_type:
        pmd_top_atomtypes = _atom_types_from_pmd(structure)
        # Consolidate parmed bondtypes and relate to topology bondtypes
        pmd_top_bondtypes = _bond_types_from_pmd(structure)
        # Consolidate parmed angletypes and relate to topology angletypes
        pmd_top_angletypes = _angle_types_from_pmd(structure)
        # Consolidate parmed dihedraltypes and relate to topology dihedraltypes
        pmd_top_dihedraltypes = _dihedral_types_from_pmd(structure)
        pmd_top_rbtorsiontypes = _rbtorsion_types_from_pmd(structure)
        # Consolidate parmed dihedral/impropertypes and relate to topology impropertypes
        pmd_top_periodic_impropertypes = _improper_types_periodic_from_pmd(
            structure
        )
        pmd_top_harmonic_impropertypes = _improper_types_harmonic_from_pmd(
            structure
        )

    ind_res = _check_independent_residues(structure)
    for residue in structure.residues:
        for atom in residue.atoms:
            element = (
                element_by_atomic_number(atom.element) if atom.element else None
            )
            site = gmso.Atom(
                name=atom.name,
                charge=atom.charge * u.elementary_charge,
                position=[atom.xx, atom.xy, atom.xz] * u.angstrom,
                atom_type=None,
                residue=(residue.name, residue.idx),
                element=element,
            )
            site.molecule = (residue.name, residue.idx) if ind_res else None
            site.atom_type = (
                pmd_top_atomtypes[atom.atom_type]
                if refer_type and isinstance(atom.atom_type, pmd.AtomType)
                else None
            )

            site_map[atom] = site
            top.add_site(site)

    for bond in structure.bonds:
        # Generate bond parameters for BondType that gets passed
        # to Bond
        top_connection = gmso.Bond(
            connection_members=[site_map[bond.atom1], site_map[bond.atom2]]
        )
        if refer_type and isinstance(bond.type, pmd.BondType):
            key = (
                bond.type.k,
                bond.type.req,
                tuple(sorted((bond.atom1.type, bond.atom2.type))),
            )
            top_connection.bond_type = pmd_top_bondtypes[key]
        top.add_connection(top_connection, update_types=False)

    for angle in structure.angles:
        # Generate angle parameters for AngleType that gets passed
        # to Angle
        top_connection = gmso.Angle(
            connection_members=[
                site_map[angle.atom1],
                site_map[angle.atom2],
                site_map[angle.atom3],
            ]
        )
        if refer_type and isinstance(angle.type, pmd.AngleType):
            member_types = (
                angle.atom2.type,
                *sorted((angle.atom1.type, angle.atom3.type)),
            )
            key = (angle.type.k, angle.type.theteq, member_types)
            top_connection.angle_type = pmd_top_angletypes[key]
        top.add_connection(top_connection, update_types=False)

    for dihedral in structure.dihedrals:
        # Generate parameters for ImproperType or DihedralType that gets passed
        # to corresponding Dihedral or Improper
        # These all follow periodic torsions functions
        # Which are the default expression in top.DihedralType
        # These periodic torsion dihedrals get stored in top.dihedrals
        # and periodic torsion impropers get stored in top.impropers

        if dihedral.improper:
            warnings.warn(
                "ParmEd improper dihedral {} ".format(dihedral)
                + "following periodic torsion "
                + "expression detected, currently accounted for as "
                + "topology.Improper with a periodic improper expression"
            )
            # TODO: Improper atom order is not always clear in a Parmed object.
            # This reader assumes the order of impropers is central atom first,
            # so that is where the central atom is located. This decision comes
            # from .top files in utils/files/NN-dimethylformamide.top, which
            # clearly places the periodic impropers with central atom listed first,
            # and that is where the atom is placed in the parmed.dihedrals object.
            top_connection = gmso.Improper(
                connection_members=[
                    site_map[dihedral.atom1],
                    site_map[dihedral.atom2],
                    site_map[dihedral.atom3],
                    site_map[dihedral.atom4],
                ],
            )
            if refer_type and isinstance(dihedral.type, pmd.DihedralType):
                key = (
                    dihedral.type.phi_k,
                    dihedral.type.phase,
                    dihedral.type.per,
                    (atomtype_namegetter(dihedral)),
                )
                top_connection.improper_type = pmd_top_periodic_impropertypes[
                    key
                ]
        else:
            top_connection = gmso.Dihedral(
                connection_members=[
                    site_map[dihedral.atom1],
                    site_map[dihedral.atom2],
                    site_map[dihedral.atom3],
                    site_map[dihedral.atom4],
                ]
            )
            if refer_type and isinstance(dihedral.type, pmd.DihedralType):
                key = (
                    dihedral.type.phi_k,
                    dihedral.type.phase,
                    dihedral.type.per,
                    (
                        tuple(reversed(atomtype_namegetter(dihedral)))
                        if rev_dih_order(dihedral)
                        else atomtype_namegetter(dihedral)
                    ),
                )
                top_connection.dihedral_type = pmd_top_dihedraltypes[key]
            # No bond parameters, make Connection with no connection_type
        top.add_connection(top_connection, update_types=False)

    for rb_torsion in structure.rb_torsions:
        # Generate dihedral parameters for DihedralType that gets passed
        # to Dihedral
        # These all follow RB torsion functions
        # These RB torsion dihedrals get stored in top.dihedrals
        if rb_torsion.improper:
            warnings.warn(
                "ParmEd improper dihedral {} ".format(rb_torsion)
                + "following RB torsion "
                + "expression detected, currently accounted for as "
                + "topology.Dihedral with a RB torsion expression"
            )

        top_connection = gmso.Dihedral(
            connection_members=[
                site_map[rb_torsion.atom1],
                site_map[rb_torsion.atom2],
                site_map[rb_torsion.atom3],
                site_map[rb_torsion.atom4],
            ],
        )
        if refer_type and isinstance(rb_torsion.type, pmd.RBTorsionType):
            key = (
                rb_torsion.type.c0,
                rb_torsion.type.c1,
                rb_torsion.type.c2,
                rb_torsion.type.c3,
                rb_torsion.type.c4,
                rb_torsion.type.c5,
                (
                    tuple(reversed(atomtype_namegetter(rb_torsion)))
                    if rev_dih_order(rb_torsion)
                    else atomtype_namegetter(rb_torsion)
                ),
            )
            top_connection.dihedral_type = pmd_top_rbtorsiontypes[key]
        top.add_connection(top_connection, update_types=False)

    for improper in structure.impropers:
        # TODO: Improper atom order is not always clear in a Parmed object.
        # This reader assumes the order of impropers is central atom first,
        # so that is where the central atom is located. This decision comes
        # from .top files in utils/files/NN-dimethylformamide.top, which
        # clearly places the periodic impropers with central atom listed first,
        # and that is where the atom is placed in the parmed.dihedrals object.
        top_connection = gmso.Improper(
            connection_members=[
                site_map[improper.atom1],
                site_map[improper.atom2],
                site_map[improper.atom3],
                site_map[improper.atom4],
            ],
        )
        if refer_type and isinstance(improper.type, pmd.ImproperType):
            # Impropers have no sorting for orders
            key = (
                improper.type.psi_k,
                improper.type.psi_eq,
                (atomtype_namegetter(improper)),
            )
            top_connection.improper_type = pmd_top_harmonic_impropertypes[key]
        top.add_connection(top_connection, update_types=False)

    top.update_topology()
    top.combining_rule = structure.combining_rule
    return top


def _atom_types_from_pmd(structure):
    """Convert ParmEd atomtypes to GMSO AtomType.

    This function take in a Parmed Structure, iterate through its
    atom's atom_type, create a corresponding GMSO.AtomType, and
    finally return a dictionary containing all pairs of pmd.AtomType
    and GMSO.AtomType

    Parameter
    ----------
        structure: pmd.Structure
            Parmed Structure that needed to be converted.

    Return
    ------
        pmd_top_atomtypes : dict
            A dictionary linking a pmd.AtomType object to its
            corresponding GMSO.AtomType object.
    """
    unique_atom_types = set()
    for atom in structure.atoms:
        if isinstance(atom.atom_type, pmd.AtomType):
            unique_atom_types.add(atom.atom_type)
    unique_atom_types = list(unique_atom_types)
    pmd_top_atomtypes = {}
    for atom_type in unique_atom_types:
        if atom_type.atomic_number:
            element = element_by_atomic_number(atom_type.atomic_number).symbol
        else:
            element = atom_type.name
        top_atomtype = gmso.AtomType(
            name=atom_type.name,
            charge=atom_type.charge * u.elementary_charge,
            tags={"element": element},
            expression="4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            parameters={
                "sigma": atom_type.sigma * u.angstrom,
                "epsilon": atom_type.epsilon * u.Unit("kcal / mol"),
            },
            independent_variables={"r"},
            mass=atom_type.mass,
        )
        pmd_top_atomtypes[atom_type] = top_atomtype
    return pmd_top_atomtypes


def _bond_types_from_pmd(structure):
    """Convert ParmEd bondtypes to GMSO BondType.

    This function takes in a Parmed Structure, iterate through its
    bond_types, create a corresponding GMSO.BondType, and finally
    return a dictionary containing all pairs of pmd.BondType
    and GMSO.BondType
    Parameters
    ----------
    structure: pmd.Structure
        Parmed Structure that needed to be converted.

    Returns
    -------
    pmd_top_bondtypes : dict
        A dictionary linking a pmd.BondType object to its
        corresponding GMSO.BondType object.
    """
    pmd_top_bondtypes = dict()
    harmonicbond_potential = lib["HarmonicBondPotential"]
    name = harmonicbond_potential.name
    expression = harmonicbond_potential.expression
    variables = harmonicbond_potential.independent_variables

    if not structure.bond_types:
        return pmd_top_bondtypes
    unique_bond_types = list(
        dict.fromkeys(
            [
                (
                    bond.type.k,
                    bond.type.req,
                    tuple(sorted((bond.atom1.type, bond.atom2.type))),
                )
                for bond in structure.bonds
            ]
        )
    )
    for btype in unique_bond_types:
        bond_params = {
            "k": (2 * btype[0] * u.Unit("kcal / (angstrom**2 * mol)")),
            "r_eq": btype[1] * u.angstrom,
        }
        member_types = btype[2]
        top_bondtype = gmso.BondType(
            name=name,
            parameters=bond_params,
            expression=expression,
            independent_variables=variables,
            member_types=member_types,
        )
        pmd_top_bondtypes[btype] = top_bondtype
    return pmd_top_bondtypes


def _angle_types_from_pmd(structure):
    """Convert ParmEd angle types to  GMSO AngleType.

    This function takes in a Parmed Structure, iterates through its
    angle_types, create a corresponding GMSO.AngleType, and finally
    return a dictionary containing all pairs of pmd.AngleType
    and GMSO.AngleType
    Parameters
    ----------
    structure: pmd.Structure
        Parmed Structure that needed to be converted.

    Returns
    -------
    pmd_top_angletypes : dict
        A dictionary linking a pmd.AngleType object to its
        corresponding GMSO.AngleType object.
    """
    pmd_top_angletypes = dict()
    harmonicbond_potential = lib["HarmonicAnglePotential"]
    name = harmonicbond_potential.name
    expression = harmonicbond_potential.expression
    variables = harmonicbond_potential.independent_variables

    if not structure.angle_types:
        return pmd_top_angletypes
    unique_angle_types = list(
        dict.fromkeys(
            [
                (
                    angle.type.k,
                    angle.type.theteq,
                    (
                        angle.atom2.type,
                        *sorted((angle.atom1.type, angle.atom3.type)),
                    ),
                )
                for angle in structure.angles
            ]
        )
    )
    for angletype in unique_angle_types:
        angle_params = {
            "k": (2 * angletype[0] * u.Unit("kcal / (rad**2 * mol)")),
            "theta_eq": (angletype[1] * u.degree),
        }
        # TODO: we need to worry about Urey Bradley terms
        # For Urey Bradley:
        # k in (kcal/(angstrom**2 * mol))
        # r_eq in angstrom
        top_angletype = gmso.AngleType(
            name=name,
            parameters=angle_params,
            expression=expression,
            independent_variables=variables,
            member_types=angletype[2],
        )
        pmd_top_angletypes[angletype] = top_angletype
    return pmd_top_angletypes


def _dihedral_types_from_pmd(structure):
    """Convert ParmEd dihedral types to GMSO DihedralType.

    This function take in a Parmed Structure, iterate through its
    dihedral_types, create a corresponding
    GMSO.DihedralType, and finally return a dictionary containing all
    pairs of pmd.Dihedraltype and GMSO.DihedralType
    Parameters
    ----------
    structure: pmd.Structure
        Parmed Structure that needed to be converted.

    Returns
    -------
    pmd_top_dihedraltypes : dict
        A dictionary linking a pmd.DihedralType
        object to its corresponding GMSO.DihedralType object.
    """
    pmd_top_dihedraltypes = dict()
    proper_dihedralsList = [
        dihedral for dihedral in structure.dihedrals if not dihedral.improper
    ]
    if len(proper_dihedralsList) == 0 or not structure.dihedral_types:
        return pmd_top_dihedraltypes
    unique_dihedral_types = list(
        dict.fromkeys(
            [
                (
                    dihedral.type.phi_k,
                    dihedral.type.phase,
                    dihedral.type.per,
                    (
                        tuple(reversed(atomtype_namegetter(dihedral)))
                        if rev_dih_order(dihedral)
                        else atomtype_namegetter(dihedral)
                    ),
                )
                for dihedral in proper_dihedralsList
            ]
        )
    )
    for dihedraltype in unique_dihedral_types:
        dihedral_params = {
            "k": (dihedraltype[0] * u.Unit("kcal / mol")),
            "phi_eq": (dihedraltype[1] * u.degree),
            "n": dihedraltype[2] * u.dimensionless,
        }
        periodic_torsion_potential = lib["PeriodicTorsionPotential"]
        name = periodic_torsion_potential.name
        expression = periodic_torsion_potential.expression
        variables = periodic_torsion_potential.independent_variables

        top_dihedraltype = gmso.DihedralType(
            name=name,
            parameters=dihedral_params,
            expression=expression,
            independent_variables=variables,
            member_types=dihedraltype[3],
        )
        pmd_top_dihedraltypes[dihedraltype] = top_dihedraltype

    return pmd_top_dihedraltypes


def _rbtorsion_types_from_pmd(structure):
    """Convert ParmEd rb_torsion types to GMSO DihedralType.

    This function take in a Parmed Structure, iterate through its
    rb_torsion_types and rb_torsion_types, create a corresponding
    GMSO.DihedralType, and finally return a dictionary containing all
    pairs of pmd.RBTorsionType and GMSO.DihedralType
    Parameters
    ----------
    structure: pmd.Structure
        Parmed Structure that needed to be converted.

    Returns
    -------
    pmd_top_dihedraltypes : dict
        A dictionary linking a pmd.RBTorsionType
        object to its corresponding GMSO.DihedralType object.
    """
    pmd_top_rbtorsiontypes = dict()
    if not structure.rb_torsion_types:
        return pmd_top_rbtorsiontypes
    unique_rb_types = list(
        dict.fromkeys(
            [
                (
                    dihedral.type.c0,
                    dihedral.type.c1,
                    dihedral.type.c2,
                    dihedral.type.c3,
                    dihedral.type.c4,
                    dihedral.type.c5,
                    (
                        tuple(reversed(atomtype_namegetter(dihedral)))
                        if rev_dih_order(dihedral)
                        else atomtype_namegetter(dihedral)
                    ),
                )
                for dihedral in structure.rb_torsions
            ]
        )
    )

    for dihedraltype in unique_rb_types:
        dihedral_params = {
            "c0": (dihedraltype[0] * u.Unit("kcal/mol")),
            "c1": (dihedraltype[1] * u.Unit("kcal/mol")),
            "c2": (dihedraltype[2] * u.Unit("kcal/mol")),
            "c3": (dihedraltype[3] * u.Unit("kcal/mol")),
            "c4": (dihedraltype[4] * u.Unit("kcal/mol")),
            "c5": (dihedraltype[5] * u.Unit("kcal/mol")),
        }

        ryckaert_bellemans_torsion_potential = lib[
            "RyckaertBellemansTorsionPotential"
        ]
        name = ryckaert_bellemans_torsion_potential.name
        expression = ryckaert_bellemans_torsion_potential.expression
        variables = ryckaert_bellemans_torsion_potential.independent_variables

        top_dihedraltype = gmso.DihedralType(
            name=name,
            parameters=dihedral_params,
            expression=expression,
            independent_variables=variables,
            member_types=dihedraltype[6],
        )
        pmd_top_rbtorsiontypes[dihedraltype] = top_dihedraltype
    return pmd_top_rbtorsiontypes


def _improper_types_periodic_from_pmd(structure):
    """Convert ParmEd DihedralTypes to GMSO ImproperType.

    This function take in a Parmed Structure, iterate through its
    dihedral_types with the improper flag, create a corresponding
    GMSO.improperType, and finally return a dictionary containing all
    pairs of pmd.impropertype and GMSO.improperType
    Parameters
    ----------
    structure: pmd.Structure
        Parmed Structure that needed to be converted.

    Returns
    -------
    pmd_top_impropertypes : dict
        A dictionary linking a pmd.improperType
        object to its corresponding GMSO.ImproperType object.
    """
    pmd_top_impropertypes = dict()
    improper_dihedralsList = [
        dihedral for dihedral in structure.dihedrals if dihedral.improper
    ]
    if len(improper_dihedralsList) == 0 or not structure.dihedral_types:
        return pmd_top_impropertypes
    unique_improper_types = list(
        dict.fromkeys(
            [
                (
                    improper.type.phi_k,
                    improper.type.phase,
                    improper.type.per,
                    (atomtype_namegetter(improper)),
                )
                for improper in improper_dihedralsList
            ]
        )
    )
    for impropertype in unique_improper_types:
        improper_params = {
            "k": (impropertype[0] * u.Unit("kcal / mol")),
            "phi_eq": (impropertype[1] * u.degree),
            "n": impropertype[2] * u.dimensionless,
        }
        periodic_torsion_potential = lib["PeriodicImproperPotential"]
        name = periodic_torsion_potential.name
        expression = periodic_torsion_potential.expression
        variables = periodic_torsion_potential.independent_variables

        top_impropertype = gmso.ImproperType(
            name=name,
            parameters=improper_params,
            expression=expression,
            independent_variables=variables,
            member_types=impropertype[3],
        )
        pmd_top_impropertypes[impropertype] = top_impropertype

    return pmd_top_impropertypes


def _improper_types_harmonic_from_pmd(
    structure,
):
    """Convert ParmEd improper types to GMSO improperType.

    This function take in a Parmed Structure, iterate through its
    improper_types, create a corresponding
    GMSO.improperType, and finally return a dictionary containing all
    pairs of pmd.impropertype and GMSO.improperType
    Parameters
    ----------
    structure: pmd.Structure
        Parmed Structure that needed to be converted.

    Returns
    -------
    pmd_top_impropertypes : dict
        A dictionary linking a pmd.improperType
        object to its corresponding GMSO.improperType object.
    """
    pmd_top_impropertypes = dict()
    if not structure.improper_types:
        return pmd_top_impropertypes
    unique_improper_types = list(
        dict.fromkeys(
            [
                (
                    improper.type.psi_k,
                    improper.type.psi_eq,
                    (atomtype_namegetter(improper)),
                )
                for improper in structure.impropers
            ]
        )
    )
    for impropertype in unique_improper_types:
        improper_params = {
            "k": (impropertype[0] * u.kcal / (u.mol * u.radian**2)),
            "phi_eq": (impropertype[1] * u.degree),
        }
        periodic_torsion_potential = lib["HarmonicTorsionPotential"]
        name = periodic_torsion_potential.name
        expression = periodic_torsion_potential.expression
        variables = periodic_torsion_potential.independent_variables

        top_impropertype = gmso.ImproperType(
            name=name,
            parameters=improper_params,
            expression=expression,
            independent_variables=variables,
            member_types=impropertype[2],
        )
        pmd_top_impropertypes[impropertype] = top_impropertype

    return pmd_top_impropertypes


def to_parmed(top, refer_type=True):
    """Convert a gmso.topology.Topology to a parmed.Structure.

    At this point we only assume a three level structure for topology
    Topology - Molecule - Residue - Sites, which transform to three level of
    Parmed Structure - Residue - Atoms (gmso Molecule level will be skipped).

    Parameters
    ----------
    top : topology.Topology
        topology.Topology instance that need to be converted
    refer_type : bool, optional, default=True
        Whether or not to transfer AtomType, BondType, AngleTye,
        and DihedralType information

    Returns
    -------
    structure : parmed.Structure
    """
    # Sanity check
    msg = "Provided argument is not a topology.Topology."
    assert isinstance(top, gmso.Topology)

    # Set up Parmed structure and define general properties
    structure = pmd.Structure()
    structure.title = top.name
    structure.box = (
        np.concatenate(
            (
                top.box.lengths.to("angstrom").value,
                top.box.angles.to("degree").value,
            )
        )
        if top.box
        else None
    )

    # Maps
    atom_map = dict()  # Map site to atom
    bond_map = dict()  # Map top's bond to structure's bond
    angle_map = dict()  # Map top's angle to strucutre's angle
    dihedral_map = dict()  # Map top's dihedral to structure's dihedral

    # Set up unparametrized system
    # Build up atom
    for site in top.sites:
        if site.element:
            atomic_number = site.element.atomic_number
        else:
            atomic_number = 0
        pmd_atom = pmd.Atom(
            atomic_number=atomic_number,
            name=site.name,
            mass=site.mass.to(u.amu).value if site.mass else None,
            charge=site.charge.to(u.elementary_charge).value
            if site.charge
            else None,
        )
        pmd_atom.xx, pmd_atom.xy, pmd_atom.xz = site.position.to(
            "angstrom"
        ).value

        # Add atom to structure
        if site.residue:
            structure.add_atom(
                pmd_atom, resname=site.residue.name, resnum=site.residue.number
            )
        else:
            structure.add_atom(pmd_atom, resname="RES", resnum=-1)
        atom_map[site] = pmd_atom

    # "Claim" all of the item it contains and subsequently index all of its item
    structure.residues.claim()

    # Create and add bonds to Parmed structure
    for bond in top.bonds:
        site1, site2 = bond.connection_members
        pmd_bond = pmd.Bond(atom_map[site1], atom_map[site2])
        structure.bonds.append(pmd_bond)
        bond_map[bond] = pmd_bond

    # Create and add angles to Parmed structure
    for angle in top.angles:
        site1, site2, site3 = angle.connection_members
        pmd_angle = pmd.Angle(atom_map[site1], atom_map[site2], atom_map[site3])
        structure.angles.append(pmd_angle)
        angle_map[angle] = pmd_angle

    # Create and add dihedrals to Parmed structure

    for dihedral in top.dihedrals:
        site1, site2, site3, site4 = dihedral.connection_members
        pmd_dihedral = pmd.Dihedral(
            atom_map[site1], atom_map[site2], atom_map[site3], atom_map[site4]
        )
        if (
            dihedral.connection_type
            and dihedral.connection_type.expression
            == parse_expr(
                "c0 * cos(phi)**0 + "
                + "c1 * cos(phi)**1 + "
                + "c2 * cos(phi)**2 + "
                + "c3 * cos(phi)**3 + "
                + "c4 * cos(phi)**4 + "
                + "c5 * cos(phi)**5"
            )
        ):
            structure.rb_torsions.append(pmd_dihedral)
        else:
            structure.dihedrals.append(pmd_dihedral)
        dihedral_map[dihedral] = pmd_dihedral

    # Set up structure for Connection Type conversion
    if refer_type:
        # Need to add a warning if Topology does not have types information
        if top.atom_types:
            _atom_types_from_gmso(top, structure, atom_map)
        if top.bond_types:
            _bond_types_from_gmso(top, structure, bond_map)
        if top.angle_types:
            _angle_types_from_gmso(top, structure, angle_map)
        if top.dihedral_types:
            _dihedral_types_from_gmso(top, structure, dihedral_map)

    return structure


def _check_independent_residues(structure):
    """Check to see if residues will constitute independent graphs."""
    # Copy from foyer forcefield.py
    for res in structure.residues:
        atoms_in_residue = set([*res.atoms])
        bond_partners_in_residue = [
            item
            for sublist in [atom.bond_partners for atom in res.atoms]
            for item in sublist
        ]
        # Handle the case of a 'residue' with no neighbors
        if not bond_partners_in_residue:
            continue
        if set(atoms_in_residue) != set(bond_partners_in_residue):
            return False
    return True


def _atom_types_from_gmso(top, structure, atom_map):
    """Convert gmso.Topology AtomType to parmed.Structure AtomType.

    This function will first check the AtomType expression of Topology and make sure it match with the one default in Parmed.
    After that, it would start atomtyping and parametrizing this part of the structure.

    Parameters
    ----------
    top : topology.Topology
        The topology that need to be converted
    structure: parmed.Structure
        The destination parmed Structure
    """
    # Maps
    atype_map = dict()
    for atom_type in top.atom_types:
        msg = "Atom type {} expression does not match Parmed AtomType default expression".format(
            atom_type.name
        )
        assert atom_type.expression == parse_expr(
            "4*epsilon*(-sigma**6/r**6 + sigma**12/r**12)"
        ), msg
        # Extract Topology atom type information
        atype_name = atom_type.name
        # Convert charge to elementary_charge
        atype_charge = float(atom_type.charge.to("Coulomb").value) / (
            1.6 * 10 ** (-19)
        )
        atype_sigma = float(atom_type.parameters["sigma"].to("angstrom").value)
        atype_epsilon = float(
            atom_type.parameters["epsilon"].to("kcal/mol").value
        )
        atype_element = element_by_atom_type(atom_type)
        atype_rmin = atype_sigma * 2 ** (1 / 6) / 2  # to rmin/2
        # Create unique Parmed AtomType object
        atype = pmd.AtomType(
            atype_name,
            None,
            atype_element.mass,
            atype_element.atomic_number,
            atype_charge,
        )
        atype.set_lj_params(atype_epsilon, atype_rmin)
        # Type map to match AtomType to its name
        atype_map[atype_name] = atype

    for site in top.sites:
        # Assign atom_type to atom
        pmd_atom = atom_map[site]
        pmd_atom.type = site.atom_type.name
        pmd_atom.atom_type = atype_map[site.atom_type.name]


def _bond_types_from_gmso(top, structure, bond_map):
    """Convert gmso.Topology BondType to parmed.Structure BondType.

    This function will first check the BondType expression of Topology and make sure it match with the one default in Parmed.
    After that, it would start atomtyping and parametrizing this part of the structure.

    Parameters
    ----------
    top : topology.Topology
        The topology that need to be converted
    structure: parmed.Structure
        The destination parmed Structure
    """
    btype_map = dict()
    for bond_type in top.bond_types:
        msg = "Bond type {} expression does not match Parmed BondType default expression".format(
            bond_type.name
        )
        assert bond_type.expression == parse_expr("0.5 * k * (r-r_eq)**2"), msg
        # Extract Topology bond_type information
        btype_k = 0.5 * float(
            bond_type.parameters["k"].to("kcal / (angstrom**2 * mol)").value
        )
        btype_r_eq = float(bond_type.parameters["r_eq"].to("angstrom").value)
        # Create unique Parmed BondType object
        btype = pmd.BondType(btype_k, btype_r_eq)
        # Type map to match Topology BondType with Parmed BondType
        btype_map[bond_type] = btype
        # Add BondType to structure.bond_types
        structure.bond_types.append(btype)

    for bond in top.bonds:
        # Assign bond_type to bond
        pmd_bond = bond_map[bond]
        pmd_bond.type = btype_map[bond.connection_type]
    structure.bond_types.claim()


def _angle_types_from_gmso(top, structure, angle_map):
    """Convert gmso.Topology AngleType to parmed.Structure AngleType.

    This function will first check the AngleType expression of Topology and make sure it match with the one default in Parmed.
    After that, it would start atomtyping and parametrizing the structure.

    Parameters
    ----------
    top : topology.Topology
        The topology that need to be converted
    structure: parmed.Structure
        The destination parmed Structure
    """
    agltype_map = dict()
    for angle_type in top.angle_types:
        msg = "Angle type {} expression does not match Parmed AngleType default expression".format(
            angle_type.name
        )
        assert angle_type.expression == parse_expr(
            "0.5 * k * (theta-theta_eq)**2"
        ), msg
        # Extract Topology angle_type information
        agltype_k = 0.5 * float(
            angle_type.parameters["k"].to("kcal / (rad**2 * mol)").value
        )
        agltype_theta_eq = float(
            angle_type.parameters["theta_eq"].to("degree").value
        )
        # Create unique Parmed AngleType object
        agltype = pmd.AngleType(agltype_k, agltype_theta_eq)
        # Type map to match Topology AngleType with Parmed AngleType
        #
        for key, value in agltype_map.items():
            if value == agltype:
                agltype = value
                break
        agltype_map[angle_type] = agltype
        # Add AngleType to structure.angle_types
        if agltype not in structure.angle_types:
            structure.angle_types.append(agltype)

    for angle in top.angles:
        # Assign angle_type to angle
        pmd_angle = angle_map[angle]
        pmd_angle.type = agltype_map[angle.connection_type]
    structure.angle_types.claim()


def _dihedral_types_from_gmso(top, structure, dihedral_map):
    """Convert gmso.Topology DihedralType to parmed.Structure DihedralType.

    This function will first check the DihedralType expression of Topology and
    make sure it match with the one default in Parmed.
    After that, it would start atomtyping and parametrizing the structure.

    Parameters
    ----------
    top : topology.Topology
        The topology that need to be converted
    structure: parmed.Structure
        The destination parmed Structure
    """
    dtype_map = dict()
    for dihedral_type in top.dihedral_types:
        msg = "Dihedral type {} expression does not match Parmed DihedralType default expressions (Periodics, RBTorsions)".format(
            dihedral_type.name
        )
        if dihedral_type.expression == parse_expr(
            "k * (1 + cos(n * phi - phi_eq))**2"
        ):
            dtype_k = float(dihedral_type.parameters["k"].to("kcal/mol").value)
            dtype_phi_eq = float(
                dihedral_type.parameters["phi_eq"].to("degrees").value
            )
            dtype_n = float(dihedral_type.parameters["n"].value)
            # Create unique Parmed DihedralType object
            dtype = pmd.DihedralType(dtype_k, dtype_n, dtype_phi_eq)
            # Add DihedralType to structure.dihedral_types
            structure.dihedral_types.append(dtype)
        elif dihedral_type.expression == parse_expr(
            "c0 * cos(phi)**0 + "
            + "c1 * cos(phi)**1 + "
            + "c2 * cos(phi)**2 + "
            + "c3 * cos(phi)**3 + "
            + "c4 * cos(phi)**4 + "
            + "c5 * cos(phi)**5"
        ):
            dtype_c0 = float(
                dihedral_type.parameters["c0"].to("kcal/mol").value
            )
            dtype_c1 = float(
                dihedral_type.parameters["c1"].to("kcal/mol").value
            )
            dtype_c2 = float(
                dihedral_type.parameters["c2"].to("kcal/mol").value
            )
            dtype_c3 = float(
                dihedral_type.parameters["c3"].to("kcal/mol").value
            )
            dtype_c4 = float(
                dihedral_type.parameters["c4"].to("kcal/mol").value
            )
            dtype_c5 = float(
                dihedral_type.parameters["c5"].to("kcal/mol").value
            )
            # Create unique DihedralType object
            dtype = pmd.RBTorsionType(
                dtype_c0, dtype_c1, dtype_c2, dtype_c3, dtype_c4, dtype_c5
            )
            # Add RBTorsionType to structure.rb_torsion_types
            structure.rb_torsion_types.append(dtype)
        else:
            raise GMSOError("msg")
        dtype_map[dihedral_type] = dtype

    for dihedral in top.dihedrals:
        pmd_dihedral = dihedral_map[dihedral]
        pmd_dihedral.type = dtype_map[dihedral.connection_type]
    structure.dihedral_types.claim()
    structure.rb_torsions.claim()
