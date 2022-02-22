"""Module support for converting to/from ParmEd objects."""
import warnings

import numpy as np
import sympy as sym
import unyt as u
from sympy.parsing.sympy_parser import parse_expr

import gmso
from gmso.core.element import (
    element_by_atom_type,
    element_by_name,
    element_by_symbol,
)
from gmso.utils.io import has_parmed, import_

if has_parmed:
    pmd = import_("parmed")


def from_parmed(structure, refer_type=True):
    """Convert a parmed.Structure to a gmso.Topology.

    Convert a parametrized or un-parametrized parmed.Structure object to a topology.Topology.
    Specifically, this method maps Structure to Topology and Atom to Site.
    At this point, this method can only convert AtomType, BondType and AngleType.
    Conversion of DihedralType will be implement in the near future.

    Parameters
    ----------
    structure : parmed.Structure
        parmed.Structure instance that need to be converted.
    refer_type : bool, optional, default=True
        Whether or not to transfer AtomType, BondType, AngleType,
        and DihedralType information

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
    # TO DO: come up with a solution to deal with partially parametrized
    # Parmed Structure
    # Old code:
    # simple check if our pmd.Structure is fully parametrized
    # is_parametrized = True if (isinstance(structure.atoms[i].atom_type,
    #        pmd.AtomType) for i in range(len(structure.atoms))) else False

    # Consolidate parmed atomtypes and relate topology atomtypes
    if refer_type:
        pmd_top_atomtypes = _atom_types_from_pmd(structure)
        # Consolidate parmed bondtypes and relate to topology bondtypes
        bond_types_map = _get_types_map(structure, "bonds")
        pmd_top_bondtypes = _bond_types_from_pmd(
            structure, bond_types_members_map=bond_types_map
        )
        # Consolidate parmed angletypes and relate to topology angletypes
        angle_types_map = _get_types_map(structure, "angles")
        pmd_top_angletypes = _angle_types_from_pmd(
            structure, angle_types_member_map=angle_types_map
        )
        # Consolidate parmed dihedraltypes and relate to topology dihedraltypes
        dihedral_types_map = _get_types_map(structure, "dihedrals")
        dihedral_types_map.update(_get_types_map(structure, "rb_torsions"))
        pmd_top_dihedraltypes = _dihedral_types_from_pmd(
            structure, dihedral_types_member_map=dihedral_types_map
        )

    subtops = list()
    for residue in structure.residues:
        subtop_name = ("{}[{}]").format(residue.name, residue.idx)
        subtops.append(gmso.SubTopology(name=subtop_name, parent=top))
        for atom in residue.atoms:
            if refer_type and isinstance(atom.atom_type, pmd.AtomType):
                site = gmso.Atom(
                    name=atom.name,
                    charge=atom.charge * u.elementary_charge,
                    position=(
                        [atom.xx, atom.xy, atom.xz] * u.angstrom
                    ).in_units(u.nm),
                    atom_type=pmd_top_atomtypes[atom.atom_type],
                    residue_name=residue.name,
                    residue_number=residue.idx,
                )
            else:
                site = gmso.Atom(
                    name=atom.name,
                    charge=atom.charge * u.elementary_charge,
                    position=(
                        [atom.xx, atom.xy, atom.xz] * u.angstrom
                    ).in_units(u.nm),
                    atom_type=None,
                    residue_name=residue.name,
                    residue_number=residue.idx,
                )
            site_map[atom] = site
            subtops[-1].add_site(site)
        top.add_subtopology(subtops[-1])

    for bond in structure.bonds:
        # Generate bond parameters for BondType that gets passed
        # to Bond
        if refer_type and isinstance(bond.type, pmd.BondType):
            top_connection = gmso.Bond(
                connection_members=[site_map[bond.atom1], site_map[bond.atom2]],
                bond_type=pmd_top_bondtypes[bond.type],
            )

        # No bond parameters, make Connection with no connection_type
        else:
            top_connection = gmso.Bond(
                connection_members=[site_map[bond.atom1], site_map[bond.atom2]],
                bond_type=None,
            )

        top.add_connection(top_connection, update_types=False)
    top.update_topology()

    for angle in structure.angles:
        # Generate angle parameters for AngleType that gets passed
        # to Angle
        if refer_type and isinstance(angle.type, pmd.AngleType):
            top_connection = gmso.Angle(
                connection_members=[
                    site_map[angle.atom1],
                    site_map[angle.atom2],
                    site_map[angle.atom3],
                ],
                angle_type=pmd_top_angletypes[angle.type],
            )
        # No bond parameters, make Connection with no connection_type
        else:
            top_connection = gmso.Angle(
                connection_members=[
                    site_map[angle.atom1],
                    site_map[angle.atom2],
                    site_map[angle.atom3],
                ],
                angle_type=None,
            )
        top.add_connection(top_connection, update_types=False)

    for dihedral in structure.dihedrals:
        # Generate dihedral parameters for DihedralType that gets passed
        # to Dihedral
        # These all follow periodic torsions functions
        # (even if they are improper dihedrals)
        # Which are the default expression in top.DihedralType
        # These periodic torsion dihedrals get stored in top.dihedrals
        if dihedral.improper:
            warnings.warn(
                "ParmEd improper dihedral {} ".format(dihedral)
                + "following periodic torsion "
                + "expression detected, currently accounted for as "
                + "topology.Dihedral with a periodic torsion expression"
            )
        if refer_type and isinstance(dihedral.type, pmd.DihedralType):
            top_connection = gmso.Dihedral(
                connection_members=[
                    site_map[dihedral.atom1],
                    site_map[dihedral.atom2],
                    site_map[dihedral.atom3],
                    site_map[dihedral.atom4],
                ],
                dihedral_type=pmd_top_dihedraltypes[dihedral.type],
            )
        # No bond parameters, make Connection with no connection_type
        else:
            top_connection = gmso.Dihedral(
                connection_members=[
                    site_map[dihedral.atom1],
                    site_map[dihedral.atom2],
                    site_map[dihedral.atom3],
                    site_map[dihedral.atom4],
                ],
                dihedral_type=None,
            )
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
        if refer_type and isinstance(rb_torsion.type, pmd.RBTorsionType):
            top_connection = gmso.Dihedral(
                connection_members=[
                    site_map[rb_torsion.atom1],
                    site_map[rb_torsion.atom2],
                    site_map[rb_torsion.atom3],
                    site_map[rb_torsion.atom4],
                ],
                dihedral_type=pmd_top_dihedraltypes[rb_torsion.type],
            )
        # No bond parameters, make Connection with no connection_type
        else:
            top_connection = gmso.Dihedral(
                connection_members=[
                    site_map[rb_torsion.atom1],
                    site_map[rb_torsion.atom2],
                    site_map[rb_torsion.atom3],
                    site_map[rb_torsion.atom4],
                ],
                dihedral_type=None,
            )
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
        top_atomtype = gmso.AtomType(
            name=atom_type.name,
            charge=atom_type.charge * u.elementary_charge,
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


def _bond_types_from_pmd(structure, bond_types_members_map=None):
    """Convert ParmEd bondtypes to GMSO BondType.

    This function takes in a Parmed Structure, iterate through its
    bond_types, create a corresponding GMSO.BondType, and finally
    return a dictionary containing all pairs of pmd.BondType
    and GMSO.BondType

    Parameters
    ----------
    structure: pmd.Structure
        Parmed Structure that needed to be converted.
    bond_types_members_map: optional, dict, default=None
        The member types (atomtype string) for each atom associated with the bond_types the structure

    Returns
    -------
    pmd_top_bondtypes : dict
        A dictionary linking a pmd.BondType object to its
        corresponding GMSO.BondType object.
    """
    pmd_top_bondtypes = dict()
    bond_types_members_map = _assert_dict(
        bond_types_members_map, "bond_types_members_map"
    )
    for btype in structure.bond_types:
        bond_params = {
            "k": (2 * btype.k * u.Unit("kcal / (angstrom**2 * mol)")),
            "r_eq": btype.req * u.angstrom,
        }
        expr = gmso.BondType._default_potential_expr()
        expr.set(parameters=bond_params)

        member_types = bond_types_members_map.get(id(btype))
        top_bondtype = gmso.BondType(
            potential_expression=expr, member_types=member_types
        )
        pmd_top_bondtypes[btype] = top_bondtype
    return pmd_top_bondtypes


def _angle_types_from_pmd(structure, angle_types_member_map=None):
    """Convert ParmEd angle types to  GMSO AngleType.

    This function takes in a Parmed Structure, iterates through its
    angle_types, create a corresponding GMSO.AngleType, and finally
    return a dictionary containing all pairs of pmd.AngleType
    and GMSO.AngleType

    Parameters
    ----------
    structure: pmd.Structure
        Parmed Structure that needed to be converted.
    angle_types_member_map: optional, dict, default=None
        The member types (atomtype string) for each atom associated with the angle_types the structure

    Returns
    -------
    pmd_top_angletypes : dict
        A dictionary linking a pmd.AngleType object to its
        corresponding GMSO.AngleType object.
    """
    pmd_top_angletypes = dict()
    angle_types_member_map = _assert_dict(
        angle_types_member_map, "angle_types_member_map"
    )

    for angletype in structure.angle_types:
        angle_params = {
            "k": (2 * angletype.k * u.Unit("kcal / (rad**2 * mol)")),
            "theta_eq": (angletype.theteq * u.degree),
        }
        expr = gmso.AngleType._default_potential_expr()
        expr.parameters = angle_params
        # Do we need to worry about Urey Bradley terms
        # For Urey Bradley:
        # k in (kcal/(angstrom**2 * mol))
        # r_eq in angstrom
        member_types = angle_types_member_map.get(id(angletype))
        top_angletype = gmso.AngleType(
            potential_expression=expr, member_types=member_types
        )
        pmd_top_angletypes[angletype] = top_angletype
    return pmd_top_angletypes


def _dihedral_types_from_pmd(structure, dihedral_types_member_map=None):
    """Convert ParmEd dihedral types to GMSO DihedralType.

    This function take in a Parmed Structure, iterate through its
    dihedral_types and rb_torsion_types, create a corresponding
    GMSO.DihedralType, and finally return a dictionary containing all
    pairs of pmd.Dihedraltype (or pmd.RBTorsionType) and GMSO.DihedralType

    Parameters
    ----------
    structure: pmd.Structure
        Parmed Structure that needed to be converted.
    dihedral_types_member_map: optional, dict, default=None
        The member types (atomtype string) for each atom associated with the dihedral_types the structure

    Returns
    -------
    pmd_top_dihedraltypes : dict
        A dictionary linking a pmd.DihedralType or pmd.RBTorsionType
        object to its corresponding GMSO.DihedralType object.
    """
    pmd_top_dihedraltypes = dict()
    dihedral_types_member_map = _assert_dict(
        dihedral_types_member_map, "dihedral_types_member_map"
    )

    for dihedraltype in structure.dihedral_types:
        dihedral_params = {
            "k": (dihedraltype.phi_k * u.Unit("kcal / mol")),
            "phi_eq": (dihedraltype.phase * u.degree),
            "n": dihedraltype.per * u.dimensionless,
        }
        expr = gmso.DihedralType._default_potential_expr()
        expr.parameters = dihedral_params
        member_types = dihedral_types_member_map.get(id(dihedraltype))
        top_dihedraltype = gmso.DihedralType(
            potential_expression=expr, member_types=member_types
        )
        pmd_top_dihedraltypes[dihedraltype] = top_dihedraltype

    for dihedraltype in structure.rb_torsion_types:
        dihedral_params = {
            "c0": (dihedraltype.c0 * u.Unit("kcal/mol")),
            "c1": (dihedraltype.c1 * u.Unit("kcal/mol")),
            "c2": (dihedraltype.c2 * u.Unit("kcal/mol")),
            "c3": (dihedraltype.c3 * u.Unit("kcal/mol")),
            "c4": (dihedraltype.c4 * u.Unit("kcal/mol")),
            "c5": (dihedraltype.c5 * u.Unit("kcal/mol")),
        }

        member_types = dihedral_types_member_map.get(id(dihedraltype))

        top_dihedraltype = gmso.DihedralType(
            parameters=dihedral_params,
            expression="c0 * cos(phi)**0 + c1 * cos(phi)**1 + "
            + "c2 * cos(phi)**2 + c3 * cos(phi)**3 + c4 * cos(phi)**4 + "
            + "c5 * cos(phi)**5",
            independent_variables="phi",
            member_types=member_types,
        )
        pmd_top_dihedraltypes[dihedraltype] = top_dihedraltype
    return pmd_top_dihedraltypes


def to_parmed(top, refer_type=True):
    """Convert a gmso.topology.Topology to a parmed.Structure.

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
    structure.box = np.concatenate(
        (
            top.box.lengths.to("angstrom").value,
            top.box.angles.to("degree").value,
        )
    )

    # Maps
    subtop_map = dict()  # Map site to subtop
    atom_map = dict()  # Map site to atom
    bond_map = dict()  # Map top's bond to structure's bond
    angle_map = dict()  # Map top's angle to strucutre's angle
    dihedral_map = dict()  # Map top's dihedral to structure's dihedral

    # Set up unparametrized system
    # Build up atom
    for site in top.sites:
        if site.element:
            atomic_number = site.element.atomic_number
            charge = site.element.charge
        else:
            atomic_number = 0
            charge = 0
        pmd_atom = pmd.Atom(
            atomic_number=atomic_number,
            name=site.name,
            mass=site.mass.to(u.amu).value,
            charge=site.charge.to(u.elementary_charge).value,
        )
        pmd_atom.xx, pmd_atom.xy, pmd_atom.xz = site.position.to(
            "angstrom"
        ).value

        # Add atom to structure
        if site.residue_name:
            structure.add_atom(
                pmd_atom, resname=site.residue_name, resnum=site.residue_number
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
        agltype_map[angle_type] = agltype
        # Add AngleType to structure.angle_types
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
            raise GMSOException("msg")
        dtype_map[dihedral_type] = dtype

    for dihedral in top.dihedrals:
        pmd_dihedral = dihedral_map[dihedral]
        pmd_dihedral.type = dtype_map[dihedral.connection_type]
    structure.dihedral_types.claim()
    structure.rb_torsions.claim()


def _get_types_map(structure, attr):
    """Build `member_types` map for atoms, bonds, angles and dihedrals."""
    assert attr in {"atoms", "bonds", "angles", "dihedrals", "rb_torsions"}
    type_map = {}
    for member in getattr(structure, attr):
        conn_type_id, member_types = _get_member_types_map_for(member)
        if conn_type_id not in type_map and all(member_types):
            type_map[conn_type_id] = member_types
    return type_map


def _get_member_types_map_for(member):
    if isinstance(member, pmd.Atom):
        return id(member.atom_type), member.type
    if isinstance(member, pmd.Bond):
        return id(member.type), (member.atom1.type, member.atom2.type)
    if isinstance(member, pmd.Angle):
        return id(member.type), (
            member.atom1.type,
            member.atom2.type,
            member.atom3.type,
        )
    if isinstance(member, pmd.Dihedral):
        return id(member.type), (
            member.atom1.type,
            member.atom2.type,
            member.atom3.type,
            member.atom4.type,
        )


def _assert_dict(input_dict, param):
    """Provide default value for a dictionary and do a type check for a parameter."""
    input_dict = {} if input_dict is None else input_dict

    if not isinstance(input_dict, dict):
        raise TypeError(
            f"Expected `{param}` to be a dictionary. "
            f"Got {type(input_dict)} instead."
        )

    return input_dict
