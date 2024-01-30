"""Module support for converting to/from ParmEd objects."""

import copy
import warnings
from operator import attrgetter, itemgetter

import numpy as np
import unyt as u
from symengine import expand

import gmso
from gmso.core.element import element_by_atomic_number, element_by_symbol
from gmso.core.views import PotentialFilters, get_parameters

pfilter = PotentialFilters.UNIQUE_PARAMETERS
from gmso.exceptions import GMSOError
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.utils.io import has_parmed, import_

if has_parmed:
    pmd = import_("parmed")

lib = PotentialTemplateLibrary()


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
        # add gmso box from structure
        top.box = gmso.Box(
            (structure.box[0:3] * u.angstrom).in_units(u.nm),
            angles=u.degree * structure.box[3:6],
        )
    top.combining_rule = structure.combining_rule

    # Consolidate parmed atomtypes and relate topology atomtypes
    if refer_type:
        pmd_top_atomtypes = _atom_types_from_pmd(structure)

    ind_res = _check_independent_residues(structure)
    for residue in structure.residues:
        for atom in residue.atoms:
            # add atom to sites in gmso
            element = (
                element_by_atomic_number(atom.element) if atom.element else None
            )
            if residue.number == -1:  # use default value of 0 in GMSO
                residue_number = 0
            else:
                residue_number = residue.number - 1
            site = gmso.Atom(
                name=atom.name,
                charge=atom.charge * u.elementary_charge,
                position=[atom.xx, atom.xy, atom.xz] * u.angstrom,
                atom_type=None,
                residue=(residue.name, residue_number),
                element=element,
            )
            site.molecule = (residue.name, residue_number) if ind_res else None
            site.atom_type = (
                copy.deepcopy(pmd_top_atomtypes[atom.atom_type])
                if refer_type and isinstance(atom.atom_type, pmd.AtomType)
                else None
            )
            if site.atom_type:
                site.atom_type.charge = atom.charge * u.elementary_charge
            site_map[atom] = site
            top.add_site(site)

    harmonicbond_potential = lib["HarmonicBondPotential"]
    name = harmonicbond_potential.name
    expression = harmonicbond_potential.expression
    variables = harmonicbond_potential.independent_variables
    for bond in structure.bonds:
        # Generate bonds and harmonic parameters
        # If typed, assumed to be harmonic bonds
        top_connection = gmso.Bond(
            connection_members=_sort_bond_members(
                top, site_map, *attrgetter("atom1", "atom2")(bond)
            )
        )
        if refer_type and isinstance(bond.type, pmd.BondType):
            conn_params = {
                "k": (2 * bond.type.k * u.Unit("kcal / (angstrom**2 * mol)")),
                "r_eq": bond.type.req * u.angstrom,
            }
            _add_conn_type_from_pmd(
                connStr="BondType",
                pmd_conn=bond,
                gmso_conn=top_connection,
                conn_params=conn_params,
                name=name,
                expression=expression,
                variables=variables,
            )
        top.add_connection(top_connection, update_types=False)

    harmonicangle_potential = lib["HarmonicAnglePotential"]
    name = harmonicangle_potential.name
    expression = harmonicangle_potential.expression
    variables = harmonicangle_potential.independent_variables
    for angle in structure.angles:
        # Generate angles and harmonic parameters
        # If typed, assumed to be harmonic angles
        top_connection = gmso.Angle(
            connection_members=_sort_angle_members(
                top, site_map, *attrgetter("atom1", "atom2", "atom3")(angle)
            )
        )
        if refer_type and isinstance(angle.type, pmd.AngleType):
            conn_params = {
                "k": (2 * angle.type.k * u.Unit("kcal / (radian**2 * mol)")),
                "theta_eq": (angle.type.theteq * u.degree),
            }
            _add_conn_type_from_pmd(
                connStr="AngleType",
                pmd_conn=angle,
                gmso_conn=top_connection,
                conn_params=conn_params,
                name=name,
                expression=expression,
                variables=variables,
            )
        top.add_connection(top_connection, update_types=False)

    periodic_torsion_potential = lib["PeriodicTorsionPotential"]
    name_proper = periodic_torsion_potential.name
    expression_proper = periodic_torsion_potential.expression
    variables_proper = periodic_torsion_potential.independent_variables
    periodic_imp_potential = lib["PeriodicImproperPotential"]
    name_improper = periodic_imp_potential.name
    expression_improper = periodic_imp_potential.expression
    variables_improper = periodic_imp_potential.independent_variables
    for dihedral in structure.dihedrals:
        # Generate dihedrals and impropers from structure.dihedrals
        # If typed, assumed to be periodic
        if dihedral.improper:
            top_connection = gmso.Improper(
                connection_members=_sort_improper_members(
                    top,
                    site_map,
                    *attrgetter("atom1", "atom2", "atom3", "atom4")(dihedral),
                )
            )
            if refer_type and isinstance(dihedral.type, pmd.DihedralType):
                conn_params = {
                    "k": (dihedral.type.phi_k * u.Unit("kcal / mol")),
                    "phi_eq": (dihedral.type.phase * u.degree),
                    "n": dihedral.type.per * u.dimensionless,
                }
                _add_conn_type_from_pmd(
                    connStr="ImproperType",
                    pmd_conn=dihedral,
                    gmso_conn=top_connection,
                    conn_params=conn_params,
                    name=name_improper,
                    expression=expression_improper,
                    variables=variables_improper,
                )
        else:
            top_connection = gmso.Dihedral(
                connection_members=_sort_dihedral_members(
                    top,
                    site_map,
                    *attrgetter("atom1", "atom2", "atom3", "atom4")(dihedral),
                )
            )
            if refer_type and isinstance(dihedral.type, pmd.DihedralType):
                conn_params = {
                    "k": (dihedral.type.phi_k * u.Unit("kcal / mol")),
                    "phi_eq": (dihedral.type.phase * u.degree),
                    "n": dihedral.type.per * u.dimensionless,
                }
                _add_conn_type_from_pmd(
                    connStr="DihedralType",
                    pmd_conn=dihedral,
                    gmso_conn=top_connection,
                    conn_params=conn_params,
                    name=name_proper,
                    expression=expression_proper,
                    variables=variables_proper,
                )
        top.add_connection(top_connection, update_types=False)

    ryckaert_bellemans_torsion_potential = lib[
        "RyckaertBellemansTorsionPotential"
    ]
    name = ryckaert_bellemans_torsion_potential.name
    expression = ryckaert_bellemans_torsion_potential.expression
    variables = ryckaert_bellemans_torsion_potential.independent_variables
    for rb_torsion in structure.rb_torsions:
        # Generate dihedrals from structure rb_torsions
        # If typed, assumed to be ryckaert bellemans torsions
        top_connection = gmso.Dihedral(
            connection_members=_sort_dihedral_members(
                top,
                site_map,
                *attrgetter("atom1", "atom2", "atom3", "atom4")(rb_torsion),
            )
        )
        if refer_type and isinstance(rb_torsion.type, pmd.RBTorsionType):
            conn_params = {
                "c0": (rb_torsion.type.c0 * u.Unit("kcal/mol")),
                "c1": (rb_torsion.type.c1 * u.Unit("kcal/mol")),
                "c2": (rb_torsion.type.c2 * u.Unit("kcal/mol")),
                "c3": (rb_torsion.type.c3 * u.Unit("kcal/mol")),
                "c4": (rb_torsion.type.c4 * u.Unit("kcal/mol")),
                "c5": (rb_torsion.type.c5 * u.Unit("kcal/mol")),
            }
            _add_conn_type_from_pmd(
                connStr="DihedralType",
                pmd_conn=rb_torsion,
                gmso_conn=top_connection,
                conn_params=conn_params,
                name=name,
                expression=expression,
                variables=variables,
            )
        top.add_connection(top_connection, update_types=False)

    periodic_torsion_potential = lib["HarmonicTorsionPotential"]
    name = periodic_torsion_potential.name
    expression = periodic_torsion_potential.expression
    variables = periodic_torsion_potential.independent_variables
    for improper in structure.impropers:
        # Generate impropers from structure impropers
        # If typed, assumed to be harmonic torsions
        top_connection = gmso.Improper(
            connection_members=_sort_improper_members(
                top,
                site_map,
                *attrgetter("atom3", "atom2", "atom1", "atom4")(improper),
            )
        )
        if refer_type and isinstance(improper.type, pmd.ImproperType):
            conn_params = {
                "k": (improper.type.psi_k * u.kcal / (u.mol * u.radian**2)),
                "phi_eq": (improper.type.psi_eq * u.degree),
            }
            _add_conn_type_from_pmd(
                connStr="ImproperType",
                pmd_conn=improper,
                gmso_conn=top_connection,
                conn_params=conn_params,
                name=name,
                expression=expression,
                variables=variables,
            )
        top.add_connection(top_connection, update_types=False)

    top.update_topology()
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
    unique_atom_types = [
        atom.atom_type
        for atom in structure.atoms
        if isinstance(atom.atom_type, pmd.AtomType)
    ]
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


def _sort_bond_members(top, site_map, atom1, atom2):
    return sorted(
        [site_map[atom1], site_map[atom2]], key=lambda x: top.get_index(x)
    )


def _sort_angle_members(top, site_map, atom1, atom2, atom3):
    sorted_angles = sorted(
        [site_map[atom1], site_map[atom3]], key=lambda x: top.get_index(x)
    )
    return (sorted_angles[0], site_map[atom2], sorted_angles[1])


# function to check reversibility of dihedral type
rev_dih_order = lambda top, site_map, x, y: top.get_index(
    site_map[x]
) > top.get_index(site_map[y])


def _sort_dihedral_members(top, site_map, atom1, atom2, atom3, atom4):
    if rev_dih_order(top, site_map, atom2, atom3):
        return itemgetter(atom4, atom3, atom2, atom1)(site_map)
    return itemgetter(atom1, atom2, atom3, atom4)(site_map)


def _sort_improper_members(top, site_map, atom1, atom2, atom3, atom4):
    sorted_impropers = sorted(
        [site_map[atom2], site_map[atom3], site_map[atom4]],
        key=lambda x: top.get_index(x),
    )
    return (site_map[atom1], *sorted_impropers)


def _add_conn_type_from_pmd(
    connStr, pmd_conn, gmso_conn, conn_params, name, expression, variables
):
    """Create a GMSO connection type and add to the conneciton object.

    This function creates the connection type object and add it to the
    connection object provided.

    Parameters
    ----------
    connStr : str
        The name of the connection type. Accepted values include
        "BondType", "AngleType", "DihedralType", and "ImproperType".
    pmd_conn : pmd.Bond/Angle/Dihedral/Improper
        The parmed connection object.
    gmso_conn : gmso.Bond/Angle/Dihedral/Improper
        The GMSO connection object.
    conn_params : dict
        The potential expression parameters in dictionary form.
    name : str
        Name of the potential form.
    expression :  expression
        The potential expression form.
    variables : dict
        The independent variables.
    """
    try:
        member_types = list(
            map(lambda x: x.atom_type.name, gmso_conn.connection_members)
        )
    except AttributeError:
        member_types = list(
            map(lambda x: f"{x}: {x.atom_type})", gmso_conn.connection_members)
        )
        raise AttributeError(
            f"Parmed structure is missing atomtypes. One of the atomtypes in \
            {member_types} is missing a type from the ParmEd structure.\
            Try using refer_type=False to not look for a parameterized structure."
        )
    try:
        get_classes = lambda x: (
            x.atom_type.atomclass if x.atom_type.atomclass else x.atom_type.name
        )
        member_classes = list(map(get_classes, gmso_conn.connection_members))
    except AttributeError:
        member_classes = list(
            map(
                lambda x: f"{x}: {x.atom_type.name})",
                gmso_conn.connection_members,
            )
        )
    top_conntype = getattr(gmso, connStr)(
        name=name,
        parameters=conn_params,
        expression=expression,
        independent_variables=variables,
        member_types=member_types,
        member_classes=member_classes,
    )
    conntypeStr = connStr.lower()[:-4] + "_type"
    setattr(gmso_conn, conntypeStr, top_conntype)


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
            charge=(
                site.charge.to(u.elementary_charge).value
                if site.charge
                else None
            ),
        )
        pmd_atom.xx, pmd_atom.xy, pmd_atom.xz = site.position.to(
            "angstrom"
        ).value

        # Add atom to structure
        if site.residue:
            structure.add_atom(
                pmd_atom,
                resname=site.residue.name,
                resnum=site.residue.number + 1,
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
        if dihedral.connection_type and expand(
            dihedral.connection_type.expression
        ) == expand(
            "c0 * cos(phi)**0 + "
            + "c1 * cos(phi)**1 + "
            + "c2 * cos(phi)**2 + "
            + "c3 * cos(phi)**3 + "
            + "c4 * cos(phi)**4 + "
            + "c5 * cos(phi)**5"
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
    for atom_type in top.atom_types(
        filter_by=PotentialFilters.UNIQUE_NAME_CLASS
    ):
        msg = "Atom type {} expression does not match Parmed AtomType default expression".format(
            atom_type.name
        )
        assert expand(atom_type.expression) == expand(
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
        if atom_type.mass:
            atype_mass = float(atom_type.mass.to("amu").value)
        else:
            atype_mass = float(
                element_by_symbol(atom_type.name).mass.to("amu").value
            )
        atype_atomic_number = getattr(
            element_by_symbol(atom_type.name), "atomic_number", None
        )
        atype_rmin = atype_sigma * 2 ** (1 / 6) / 2  # to rmin/2
        # Create unique Parmed AtomType object
        atype = pmd.AtomType(
            atype_name,
            None,
            atype_mass,
            atype_atomic_number,
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
    for bond_type in top.bond_types(filter_by=pfilter):
        msg = "Bond type {} expression does not match Parmed BondType default expression".format(
            bond_type.name
        )
        assert expand(bond_type.expression) == expand(
            "0.5 * k * (r-r_eq)**2"
        ), msg
        # Extract Topology bond_type information
        btype_k = 0.5 * float(
            bond_type.parameters["k"].to("kcal / (angstrom**2 * mol)").value
        )
        btype_r_eq = float(bond_type.parameters["r_eq"].to("angstrom").value)
        # Create unique Parmed BondType object
        btype = pmd.BondType(btype_k, btype_r_eq)
        # Type map to match Topology BondType parameters with Parmed BondType
        btype_map[get_parameters(bond_type)] = btype
        # Add BondType to structure.bond_types
        structure.bond_types.append(btype)

    for bond in top.bonds:
        # Assign bond_type to bond
        pmd_bond = bond_map[bond]
        pmd_bond.type = btype_map[get_parameters(bond.bond_type)]
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
    for angle_type in top.angle_types(filter_by=pfilter):
        msg = "Angle type {} expression does not match Parmed AngleType default expression".format(
            angle_type.name
        )
        assert expand(angle_type.expression) == expand(
            "0.5 * k * (theta-theta_eq)**2"
        ), msg
        # Extract Topology angle_type information
        agltype_k = 0.5 * float(
            angle_type.parameters["k"].to("kcal / (radian**2 * mol)").value
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
        agltype_map[get_parameters(angle_type)] = agltype
        # Add AngleType to structure.angle_types
        if agltype not in structure.angle_types:
            structure.angle_types.append(agltype)

    for angle in top.angles:
        # Assign angle_type to angle
        pmd_angle = angle_map[angle]
        pmd_angle.type = agltype_map[get_parameters(angle.connection_type)]
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
    for dihedral_type in top.dihedral_types(filter_by=pfilter):
        msg = "Dihedral type {} expression does not match Parmed DihedralType default expressions (Periodics, RBTorsions)".format(
            dihedral_type.name
        )
        if expand(dihedral_type.expression) == expand(
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
        elif expand(dihedral_type.expression) == expand(
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
                dtype_c0,
                dtype_c1,
                dtype_c2,
                dtype_c3,
                dtype_c4,
                dtype_c5,
                list=structure.rb_torsion_types,
            )
            # Add RBTorsionType to structure.rb_torsion_types
            structure.rb_torsion_types.append(dtype)
            # dtype._idx = len(structure.rb_torsion_types) - 1
        else:
            raise GMSOError("msg")
        dtype_map[get_parameters(dihedral_type)] = dtype

    for dihedral in top.dihedrals:
        pmd_dihedral = dihedral_map[dihedral]
        pmd_dihedral.type = dtype_map[get_parameters(dihedral.connection_type)]
    structure.dihedral_types.claim()
    structure.rb_torsions.claim()
