"""Write GSD files from GMSO topologies."""
from __future__ import division

import warnings

import numpy as np
import unyt as u
from unyt.array import allclose_units

from gmso.core.bond import Bond
from gmso.exceptions import NotYetImplementedWarning
from gmso.formats.formats_registry import saves_as
from gmso.utils.geometry import coord_shift
from gmso.utils.io import has_gsd

__all__ = ["write_gsd"]

if has_gsd:
    import gsd
    import gsd.hoomd


@saves_as(".gsd")
def write_gsd(
    top,
    filename,
    ref_distance=1.0 * u.nm,
    ref_mass=1.0 * u.Unit("g/mol"),
    ref_energy=1.0 * u.Unit("kcal/mol"),
    rigid_bodies=None,
    shift_coords=True,
    write_special_pairs=True,
):
    """Output a GSD file (HOOMD v2 default data format).

    The `GSD` binary file format is the native format of HOOMD-Blue. This file
    can be used as a starting point for a HOOMD-Blue simulation, for analysis,
    and for visualization in various tools.

    Parameters
    ----------
    top : gmso.Topology
        gmso.Topology object
    filename : str
        Path of the output file.
    ref_distance : float, optional, default=1.0
        Reference distance for conversion to reduced units
    ref_mass : float, optional, default=1.0
        Reference mass for conversion to reduced units
    ref_energy : float, optional, default=1.0
        Reference energy for conversion to reduced units
    rigid_bodies : list of int, optional, default=None
        List of rigid body information. An integer value is required for each
        atom corresponding to the index of the rigid body the particle is to be
        associated with. A value of None indicates the atom is not part of a
        rigid body.
    shift_coords : bool, optional, default=True
        Shift coordinates from (0, L) to (-L/2, L/2) if necessary.
    write_special_pairs : bool, optional, default=True
        Writes out special pair information necessary to correctly use the OPLS
        fudged 1,4 interactions in HOOMD.

    Notes
    -----
    Force field parameters are not written to the GSD file and must be included
    manually in a HOOMD input script. Work on a HOOMD plugin is underway to
    read force field parameters from a Foyer XML file.

    """
    xyz = u.unyt_array([site.position for site in top.sites])
    if shift_coords:
        warnings.warn("Shifting coordinates to [-L/2, L/2]")
        xyz = coord_shift(xyz, top.box)

    gsd_snapshot = gsd.hoomd.Snapshot()

    gsd_snapshot.configuration.step = 0
    gsd_snapshot.configuration.dimensions = 3

    # Write box information
    (lx, ly, lz, xy, xz, yz) = _prepare_box_information(top)
    lx = lx / ref_distance
    ly = ly / ref_distance
    lz = lz / ref_distance
    gsd_snapshot.configuration.box = np.array([lx, ly, lz, xy, xz, yz])

    warnings.warn(
        "Only writing particle and bond information."
        " Angle and dihedral is not currently written to GSD files",
        NotYetImplementedWarning,
    )
    _write_particle_information(
        gsd_snapshot, top, xyz, ref_distance, ref_mass, ref_energy, rigid_bodies
    )
    # if write_special_pairs:
    #    _write_pair_information(gsd_snapshot, top)
    if top.n_bonds > 0:
        _write_bond_information(gsd_snapshot, top)
    # if structure.angles:
    #    _write_angle_information(gsd_snapshot, top)
    # if structure.rb_torsions:
    #    _write_dihedral_information(gsd_snapshot, top)

    with gsd.hoomd.open(filename, mode="wb") as gsd_file:
        gsd_file.append(gsd_snapshot)


def _write_particle_information(
    gsd_snapshot, top, xyz, ref_distance, ref_mass, ref_energy, rigid_bodies
):
    """Write out the particle information."""
    gsd_snapshot.particles.N = top.n_sites
    warnings.warn("{} particles detected".format(top.n_sites))
    gsd_snapshot.particles.position = xyz / ref_distance

    types = [
        site.name if site.atom_type is None else site.atom_type.name
        for site in top.sites
    ]

    unique_types = list(set(types))
    unique_types = sorted(unique_types)
    gsd_snapshot.particles.types = unique_types
    warnings.warn("{} unique particle types detected".format(len(unique_types)))

    typeids = np.array([unique_types.index(t) for t in types])
    gsd_snapshot.particles.typeid = typeids

    masses = np.array([site.mass for site in top.sites])
    masses[masses == 0] = 1.0
    gsd_snapshot.particles.mass = masses / ref_mass

    charges = np.array([site.charge for site in top.sites])
    e0 = u.physical_constants.eps_0.in_units(
        u.elementary_charge**2 / u.Unit("kcal*angstrom/mol")
    )
    """
    Permittivity of free space = 2.39725e-4 e^2/((kcal/mol)(angstrom)),
    where e is the elementary charge
    """
    charge_factor = (4.0 * np.pi * e0 * ref_distance * ref_energy) ** 0.5
    gsd_snapshot.particles.charge = charges / charge_factor

    if rigid_bodies:
        warnings.warn(
            "Rigid bodies detected, but not yet implemented for GSD",
            NotYetImplementedWarning,
        )
    # if rigid_bodies:
    #    rigid_bodies = [-1 if body is None else body for body in rigid_bodies]
    # gsd_snapshot.particles.body = rigid_bodies


def _write_pair_information(gsd_snapshot, top):
    """Write the special pairs in the system.

    Parameters
    ----------
    gsd_snapshot :
        The file object of the GSD file being written
    structure : parmed.Structure
        Parmed structure object holding system information

    Warnings
    --------
    Not yet implemented for `gmso.core.topology` objects.

    """
    # pair_types = []
    # pair_typeid = []
    # pairs = []
    # for ai in structure.atoms:
    #    for aj in ai.dihedral_partners:
    #        #make sure we don't double add
    #        if ai.idx > aj.idx:
    #            ps = '-'.join(sorted([ai.type, aj.type], key=natural_sort))
    #            if ps not in pair_types:
    #                pair_types.append(ps)
    #            pair_typeid.append(pair_types.index(ps))
    #            pairs.append((ai.idx, aj.idx))
    # gsd_snapshot.pairs.types = pair_types
    # gsd_snapshot.pairs.typeid = pair_typeid
    # gsd_snapshot.pairs.group = pairs
    # gsd_snapshot.pairs.N = len(pairs)


def _write_bond_information(gsd_snapshot, top):
    """Write the bonds in the system.

    Parameters
    ----------
    gsd_snapshot :
        The file object of the GSD file being written
    top : gmso.Topology
        Topology object holding system information

    """
    gsd_snapshot.bonds.N = top.n_bonds
    warnings.warn("{} bonds detected".format(top.n_bonds))

    unique_bond_types = set()
    for bond in top.connections:
        if isinstance(bond, Bond):
            t1, t2 = (
                bond.connection_members[0].atom_type,
                bond.connection_members[1].atom_type,
            )
            if t1 is None or t2 is None:
                t1, t2 = (
                    bond.connection_members[0].name,
                    bond.connection_members[1].name,
                )
            t1, t2 = sorted([t1, t2], key=lambda x: x.name)
            bond_type = "-".join((t1.name, t2.name))

            unique_bond_types.add(bond_type)
    unique_bond_types = sorted(list(unique_bond_types))
    gsd_snapshot.bonds.types = unique_bond_types
    warnings.warn(
        "{} unique bond types detected".format(len(unique_bond_types))
    )

    bond_typeids = []
    bond_groups = []
    for bond in top.bonds:
        if isinstance(bond, Bond):
            t1, t2 = (
                bond.connection_members[0].atom_type,
                bond.connection_members[1].atom_type,
            )
            if t1 is None or t2 is None:
                t1, t2 = (
                    bond.connection_members[0].name,
                    bond.connection_members[1].name,
                )
            t1, t2 = sorted([t1, t2], key=lambda x: x.name)

            bond_type = "-".join((t1.name, t2.name))
            bond_typeids.append(unique_bond_types.index(bond_type))
            bond_groups.append(
                (
                    top.sites.index(bond.connection_members[0]),
                    top.sites.index(bond.connection_members[1]),
                )
            )

    gsd_snapshot.bonds.typeid = bond_typeids
    gsd_snapshot.bonds.group = bond_groups


def _write_angle_information(gsd_snapshot, structure):
    """Write the angles in the system.

    Parameters
    ----------
    gsd_snapshot :
        The file object of the GSD file being written
    structure : parmed.Structure
        Parmed structure object holding system information

    Warnings
    --------
    Not yet implemented for gmso.core.topology objects

    """
    # gsd_snapshot.angles.N = len(structure.angles)

    # unique_angle_types = set()
    # for angle in structure.angles:
    #    t1, t2, t3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
    #    t1, t3 = sorted([t1, t3], key=natural_sort)
    #    angle_type = ('-'.join((t1, t2, t3)))
    #    unique_angle_types.add(angle_type)
    # unique_angle_types = sorted(list(unique_angle_types), key=natural_sort)
    # gsd_snapshot.angles.types = unique_angle_types

    # angle_typeids = []
    # angle_groups = []
    # for angle in structure.angles:
    #    t1, t2, t3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
    #    t1, t3 = sorted([t1, t3], key=natural_sort)
    #    angle_type = ('-'.join((t1, t2, t3)))
    #    angle_typeids.append(unique_angle_types.index(angle_type))
    #    angle_groups.append((angle.atom1.idx, angle.atom2.idx,
    #                         angle.atom3.idx))

    # gsd_snapshot.angles.typeid = angle_typeids
    # gsd_snapshot.angles.group = angle_groups
    pass


def _write_dihedral_information(gsd_snapshot, structure):
    """Write the dihedrals in the system.

    Parameters
    ----------
    gsd_snapshot :
        The file object of the GSD file being written
    structure : parmed.Structure
        Parmed structure object holding system information

    Warnings
    --------
    Not yet implemented for gmso.core.topology objects

    """
    # gsd_snapshot.dihedrals.N = len(structure.rb_torsions)

    # unique_dihedral_types = set()
    # for dihedral in structure.rb_torsions:
    #    t1, t2 = dihedral.atom1.type, dihedral.atom2.type
    #    t3, t4 = dihedral.atom3.type, dihedral.atom4.type
    #    if [t2, t3] == sorted([t2, t3], key=natural_sort):
    #        dihedral_type = ('-'.join((t1, t2, t3, t4)))
    #    else:
    #        dihedral_type = ('-'.join((t4, t3, t2, t1)))
    #    unique_dihedral_types.add(dihedral_type)
    # unique_dihedral_types = sorted(list(unique_dihedral_types), key=natural_sort)
    # gsd_snapshot.dihedrals.types = unique_dihedral_types

    # dihedral_typeids = []
    # dihedral_groups = []
    # for dihedral in structure.rb_torsions:
    #    t1, t2 = dihedral.atom1.type, dihedral.atom2.type
    #    t3, t4 = dihedral.atom3.type, dihedral.atom4.type
    #    if [t2, t3] == sorted([t2, t3], key=natural_sort):
    #        dihedral_type = ('-'.join((t1, t2, t3, t4)))
    #    else:
    #        dihedral_type = ('-'.join((t4, t3, t2, t1)))
    #    dihedral_typeids.append(unique_dihedral_types.index(dihedral_type))
    #    dihedral_groups.append((dihedral.atom1.idx, dihedral.atom2.idx,
    #                            dihedral.atom3.idx, dihedral.atom4.idx))

    # gsd_snapshot.dihedrals.typeid = dihedral_typeids
    # gsd_snapshot.dihedrals.group = dihedral_groups
    pass


def _prepare_box_information(top):
    """Prepare the box information for writing to gsd."""
    lx = ly = lz = xy = xz = yz = 0.0
    if allclose_units(
        top.box.angles, np.array([90, 90, 90]) * u.degree, rtol=1e-5, atol=1e-8
    ):
        warnings.warn("Orthorhombic box detected")
        lx, ly, lz = top.box.lengths
        xy, xz, yz = 0.0, 0.0, 0.0
    else:
        warnings.warn("Non-orthorhombic box detected")
        u_vectors = top.box.get_unit_vectors()
        lx, ly, lz = top.box.lengths
        xy = u_vectors[1][0]
        xz = u_vectors[2][0]
        yz = u_vectors[2][1]
    return lx, ly, lz, xy, xz, yz
