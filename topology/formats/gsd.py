from __future__ import division

from collections import OrderedDict
from copy import deepcopy
from math import floor
import re

import numpy as np
import unyt as u
from oset import oset as OrderedSet
import gsd
import gsd.hoomd

from topology.core.box import Box
from topology.utils.geometry import coord_shift
from topology.testing.utils import allclose

__all__ = ['write_gsd']


def write_gsd(top,
              filename,
              ref_distance=1.0 * u.nm,
              ref_mass=1.0 * u.Unit('g/mol'),
              ref_energy=1.0 * u.Unit('kcal/mol'),
              rigid_bodies=None,
              shift_coords=True,
              write_special_pairs=True):
    """Output a GSD file (HOOMD v2 default data format).

    Parameters
    ----------
    top : topology.Topology
        topology.Topology object
    filename : str
        Path of the output file.
    ref_distance : float, optional, default=1.0
        Reference distance for conversion to reduced units
    ref_mass : float, optional, default=1.0
        Reference mass for conversion to reduced units
    ref_energy : float, optional, default=1.0
        Reference energy for conversion to reduced units
    rigid_bodies : list of int, optional, default=None
        List of rigid body information. An integer value is required for
        each atom corresponding to the index of the rigid body the particle
        is to be associated with. A value of None indicates the atom is not
        part of a rigid body.
    shift_coords : bool, optional, default=True
        Shift coordinates from (0, L) to (-L/2, L/2) if necessary.
    write_special_pairs : bool, optional, default=True
        Writes out special pair information necessary to correctly use the OPLS fudged 1,4 interactions
        in HOOMD.

    Notes
    -----
    Force field parameters are not written to the GSD file and must be included
    manually into a HOOMD input script. Work on a HOOMD plugin is underway to
    read force field parameters from a Foyer XML file.

    """

    xyz = u.unyt_array([site.position for site in top.site_list])
    if shift_coords:
        xyz = coord_shift(xyz, top.box)

    gsd_file = gsd.hoomd.Snapshot()

    gsd_file.configuration.step = 0
    gsd_file.configuration.dimensions = 3

    # Write box information
    if allclose(top.box.angles, np.array([90, 90, 90]) * u.degree):
        gsd_file.configuration.box = np.hstack((top.box.lengths / ref_distance,
                                                np.zeros(3)))
    else:
        a, b, c = top.box.lengths / ref_distance
        alpha, beta, gamma = top.box.angles

        lx = a
        xy = b * np.cos(gamma)
        xz = c * np.cos(beta)
        ly = np.sqrt(b**2 - xy**2)
        yz = (b * c * np.cos(alpha) - xy * xz) / ly
        lz = np.sqrt(c**2 - xz**2 - yz**2)

        gsd_file.configuration.box = np.array([lx, ly, lz, xy, xz, yz])

    _write_particle_information(gsd_file, top, xyz, ref_distance, ref_mass,
                                ref_energy, rigid_bodies)
    #if write_special_pairs:
    #    _write_pair_information(gsd_file, top)
    if top.n_connections > 0:
        _write_bond_information(gsd_file, top)
    #if structure.angles:
    #    _write_angle_information(gsd_file, top)
    #if structure.rb_torsions:
    #    _write_dihedral_information(gsd_file, top)

    gsd.hoomd.create(filename, gsd_file)


def _write_particle_information(gsd_file, top, xyz, ref_distance, ref_mass,
                                ref_energy, rigid_bodies):
    """Write out the particle information.

    """

    gsd_file.particles.N = top.n_sites
    gsd_file.particles.position = xyz / ref_distance

    types = [
        site.name if site.atom_type is None else site.atom_type.name
        for site in top.site_list
    ]

    unique_types = list(set(types))
    unique_types = sorted(unique_types)
    #unique_types.sort(key=natural_sort)
    gsd_file.particles.types = unique_types

    typeids = np.array([unique_types.index(t) for t in types])
    gsd_file.particles.typeid = typeids

    masses = np.array([site.mass for site in top.site_list])
    masses[masses == 0] = 1.0
    gsd_file.particles.mass = masses / ref_mass

    charges = np.array([site.charge for site in top.site_list])
    e0 = u.physical_constants.eps_0.in_units(
        u.elementary_charge**2 / u.Unit('kcal*angstrom/mol'))
    '''
    Permittivity of free space = 2.39725e-4 e^2/((kcal/mol)(angstrom)),
    where e is the elementary charge
    '''
    charge_factor = (4.0 * np.pi * e0 * ref_distance * ref_energy)**0.5
    gsd_file.particles.charge = charges / charge_factor

    if rigid_bodies:
        rigid_bodies = [-1 if body is None else body for body in rigid_bodies]
    gsd_file.particles.body = rigid_bodies


def _write_pair_information(gsd_file, top):
    """[NOT IMPLEMENTED FOR TOPOLOGY YET] Write the special pairs in the system.

        Parameters
    ----------
    gsd_file :
        The file object of the GSD file being written
    structure : parmed.Structure
        Parmed structure object holding system information
    """
    #pair_types = []
    #pair_typeid = []
    #pairs = []
    #for ai in structure.atoms:
    #    for aj in ai.dihedral_partners:
    #        #make sure we don't double add
    #        if ai.idx > aj.idx:
    #            ps = '-'.join(sorted([ai.type, aj.type], key=natural_sort))
    #            if ps not in pair_types:
    #                pair_types.append(ps)
    #            pair_typeid.append(pair_types.index(ps))
    #            pairs.append((ai.idx, aj.idx))
    #gsd_file.pairs.types = pair_types
    #gsd_file.pairs.typeid = pair_typeid
    #gsd_file.pairs.group = pairs
    #gsd_file.pairs.N = len(pairs)


def _write_bond_information(gsd_file, top):
    """Write the bonds in the system.

    Parameters
    ----------
    gsd_file :
        The file object of the GSD file being written
    top : topology.Topology
        Topology object holding system information

    """

    gsd_file.bonds.N = len(top.n_connections)

    unique_bond_types = set()
    for bond in top.connection_list:
        t1, t2 = bond.site1.atom_type, bond.site2.atom_type
        if t1 is None or t2 is None:
            t1, t2 = bond.site1.name, bond.site2.name
        t1, t2 = sorted([t1, t2])
        try:
            bond_type = ('-'.join((t1, t2)))
        except AttributeError:  # no forcefield applied, bond.type is None
            bond_type = ('-'.join((t1, t2)), 0.0, 0.0)
        unique_bond_types.add(bond_type)
    unique_bond_types = sorted(list(unique_bond_types))
    gsd_file.bonds.types = unique_bond_types

    bond_typeids = []
    bond_groups = []
    for bond in top.connections:
        t1, t2 = bond.site1.atom_type, bond.site2.atom_type
        if t1 is None or t2 is None:
            t1, t2 = bond.site1.name, bond.site2.name
        t1, t2 = sorted([t1, t2])
        try:
            bond_type = ('-'.join((t1, t2)))
        except AttributeError:  # no forcefield applied, bond.type is None
            bond_type = ('-'.join((t1, t2)), 0.0, 0.0)
        bond_typeids.append(unique_bond_types.index(bond_type))
        bond_groups.append((top.site_list.index(bond.site1),
                            top.site_list.index(bond.site2)))

    gsd_file.bonds.typeid = bond_typeids
    gsd_file.bonds.group = bond_groups


def _write_angle_information(gsd_file, structure):
    """[NOT IMPLEMENTED] Write the angles in the system.

    Parameters
    ----------
    gsd_file :
        The file object of the GSD file being written
    structure : parmed.Structure
        Parmed structure object holding system information

    """

    #gsd_file.angles.N = len(structure.angles)

    #unique_angle_types = set()
    #for angle in structure.angles:
    #    t1, t2, t3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
    #    t1, t3 = sorted([t1, t3], key=natural_sort)
    #    angle_type = ('-'.join((t1, t2, t3)))
    #    unique_angle_types.add(angle_type)
    #unique_angle_types = sorted(list(unique_angle_types), key=natural_sort)
    #gsd_file.angles.types = unique_angle_types

    #angle_typeids = []
    #angle_groups = []
    #for angle in structure.angles:
    #    t1, t2, t3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
    #    t1, t3 = sorted([t1, t3], key=natural_sort)
    #    angle_type = ('-'.join((t1, t2, t3)))
    #    angle_typeids.append(unique_angle_types.index(angle_type))
    #    angle_groups.append((angle.atom1.idx, angle.atom2.idx,
    #                         angle.atom3.idx))

    #gsd_file.angles.typeid = angle_typeids
    #gsd_file.angles.group = angle_groups


def _write_dihedral_information(gsd_file, structure):
    """[NOT IMPLEMENTED] Write the dihedrals in the system.

    Parameters
    ----------
    gsd_file :
        The file object of the GSD file being written
    structure : parmed.Structure
        Parmed structure object holding system information

    """

    #gsd_file.dihedrals.N = len(structure.rb_torsions)

    #unique_dihedral_types = set()
    #for dihedral in structure.rb_torsions:
    #    t1, t2 = dihedral.atom1.type, dihedral.atom2.type
    #    t3, t4 = dihedral.atom3.type, dihedral.atom4.type
    #    if [t2, t3] == sorted([t2, t3], key=natural_sort):
    #        dihedral_type = ('-'.join((t1, t2, t3, t4)))
    #    else:
    #        dihedral_type = ('-'.join((t4, t3, t2, t1)))
    #    unique_dihedral_types.add(dihedral_type)
    #unique_dihedral_types = sorted(list(unique_dihedral_types), key=natural_sort)
    #gsd_file.dihedrals.types = unique_dihedral_types

    #dihedral_typeids = []
    #dihedral_groups = []
    #for dihedral in structure.rb_torsions:
    #    t1, t2 = dihedral.atom1.type, dihedral.atom2.type
    #    t3, t4 = dihedral.atom3.type, dihedral.atom4.type
    #    if [t2, t3] == sorted([t2, t3], key=natural_sort):
    #        dihedral_type = ('-'.join((t1, t2, t3, t4)))
    #    else:
    #        dihedral_type = ('-'.join((t4, t3, t2, t1)))
    #    dihedral_typeids.append(unique_dihedral_types.index(dihedral_type))
    #    dihedral_groups.append((dihedral.atom1.idx, dihedral.atom2.idx,
    #                            dihedral.atom3.idx, dihedral.atom4.idx))

    #gsd_file.dihedrals.typeid = dihedral_typeids
    #gsd_file.dihedrals.group = dihedral_groups
