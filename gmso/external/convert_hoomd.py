"""Convert GMSO Topology to GSD snapshot."""

from __future__ import division

import copy
import itertools
import json
import re

import numpy as np
import unyt as u
from unyt.array import allclose_units

from gmso.core.views import PotentialFilters
from gmso.exceptions import EngineIncompatibilityError, NotYetImplementedWarning
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.utils.connectivity import generate_pairs_lists
from gmso.utils.conversions import convert_ryckaert_to_opls
from gmso.utils.geometry import coord_shift, moment_of_inertia
from gmso.utils.io import has_gsd, has_hoomd
from gmso.utils.sorting import (
    sort_by_classes,
    sort_by_types,
    sort_connection_members,
)
from gmso.utils.units import convert_params_units

if has_gsd:
    import gsd.hoomd
if has_hoomd:
    import hoomd

    hoomd_version = hoomd.version.version.split(".")
else:
    hoomd_version = None

import logging

logger = logging.getLogger(__name__)

# Note, charge will always be assumed to be in elementary_charge
MD_UNITS = {
    "energy": u.kJ / u.mol,
    "length": u.nm,
    "mass": u.g / u.mol,  # aka amu
}

AKMA_UNITS = {
    "energy": u.kcal / u.mol,
    "length": u.angstrom,
    "mass": u.g / u.mol,  # aka amu
}


def get_cell_nlist(top, buffer=0.4):
    """Create a cell neighborlist for use in hoomd-blue based on top.scaling_factors"""
    nb_scaling_factors, coul_scaling_factors = top.scaling_factors
    if all(top.scaling_factors[0] == top.scaling_factors[1]):
        outVals = None
    else:
        outVals = []  # store neighborlists
    for scalar in [nb_scaling_factors, coul_scaling_factors]:
        exclusions = list()
        for i, val in enumerate(scalar):
            if val == 1:  # skip values that are not exclusions
                continue
            if i == 0:
                exclusions.append("bond")
            else:
                exclusions.append(f"1-{i + 2}")
        if outVals is None:
            outVals = hoomd.md.nlist.Cell(exclusions=exclusions, buffer=buffer)
            return outVals, outVals  # return same neighborlist twice
        else:
            outVals.append(hoomd.md.nlist.Cell(exclusions=exclusions, buffer=buffer))
    return outVals


def to_gsd_snapshot(
    top,
    base_units=None,
    shift_coords=True,
    parse_special_pairs=True,
    auto_scale=False,
):
    """Create a gsd.hoomd.Frame objcet (HOOMD default data format).

    The gsd snapshot defines the topology of HOOMD-blue simulations. This file
    can be used as a starting point for a HOOMD-blue simulation, for analysis,
    and for visualization in various tools.

    Parameters
    ----------
    top : gmso.Topology
        gmso.Topology object
    filename : str
        Path of the output file.
    base_units : dict, optinoal, default=None
        The dictionary of base units to be converted to. Entries restricted to
        "energy", "length", and "mass". There is also option to used predefined
        unit systems ("MD" or "AKMA" provided as string). If None is provided,
        this method will perform no units conversion.
    shift_coords : bool, optional, default=True
        Shift coordinates from (0, L) to (-L/2, L/2) if necessary.
    parse_special_pairs : bool, optional, default=True
        Writes out special pair information necessary to correctly use the OPLS
        fudged 1,4 interactions in HOOMD.
    auto_scale : bool or dict, optional, default=False
        Automatically scaling relevant length, energy and mass units.
        Referenced mass unit is obtained from sites' masses.
        Referenced energy and distance are refered from sites' atom types (when applicable).
        If the referenced scaling values cannot be determined (e.g., when the topology is not typed),
        all reference scaling values is set to 1.
        A dictionary specifying the referenced scaling values may also be provided for this argument.

    Return
    ------
    gsd_snapshot : gsd.hoomd.Frame
        Converted hoomd Snapshot. Always retruned.
    base_units : dict
        Base units dictionary utilized during the conversion. Always returned.
    rigid_info : hoomd.md.constrain.Rigid
        A hoomd constraint object storing constituent particle information
        needed to run rigid body simulations in HOOMD-blue.
        This is only returned if sites with `site.molecule.isrigid = True`
        are found in the topology.

    Notes
    -----
    Force field parameters are not written to the GSD file and must be
    generated using `gmso.external.convert_hoomd.to_hoomd_forcefield()`.

    If you are using mBuild and GMSO to initialize rigid body simulations
    with HOOMD-blue set the `site.molecule.isrigid` property to `True`.
    If your topology contains a mix of rigid and flexible molecules,
    the rigid molecules must come first in the hierarchy of the GMSO Toplogy
    and therefore also in the mBuild Compound hierarchy.
    For more information about running rigid body simulations in HOOMD-blue,
    see https://hoomd-blue.readthedocs.io/en/latest/tutorial/06-Modelling-Rigid-Bodies/00-index.html

    For more information on the units from `base_units`, see
    https://hoomd-blue.readthedocs.io/en/latest/units.html#base-units

    Rigid Body Example
    ------------------
    In this example, a system with rigid benzene and flexible ethane
    is initialized. Since benzene is rigid, it must be passed first
    into `mb.fill_box` so that benzene molecules are at the top of
    `box` hierarchy.
        ::
            ethane = mb.lib.molecules.Ethane()
            ethane.name = "ethane"
            benzene = mb.load("c1ccccc1", smiles=True)
            benzene.name = "benzene"
            box = mb.fill_box([benzene, ethane], n_compounds=[1, 1], box=[2, 2, 2])
            top = from_mbuild(box)
            for site in top.sites:
                if site.molecule.name == "benzene":
                    site.molecule.isrigid = True

            snapshot, refs, rigid = to_gsd_snapshot(top)

    """
    if int(hoomd_version[0]) < 4:
        raise EngineIncompatibilityError(
            "GMSO is only compatible with HOOMD-blue >= 4.0"
        )
    base_units = _validate_base_units(base_units, top, auto_scale)
    gsd_snapshot = gsd.hoomd.Frame()

    gsd_snapshot.configuration.step = 0
    gsd_snapshot.configuration.dimensions = 3

    # Write box information
    (lx, ly, lz, xy, xz, yz) = _prepare_box_information(top)

    lx = lx.to_value(base_units["length"])
    ly = ly.to_value(base_units["length"])
    lz = lz.to_value(base_units["length"])

    gsd_snapshot.configuration.box = np.array([lx, ly, lz, xy, xz, yz])

    logger.info(
        "Only writing particle, bond, sangle, proper and improper dihedral information."
        "Special pairs are not currently written to GSD files",
    )

    # Create IndexMap
    site_indexMap = {id(site): i for i, site in enumerate(top.sites)}

    n_rigid, rigid_info = _parse_particle_information(
        gsd_snapshot,
        top,
        base_units,
        shift_coords,
        u.unyt_array([lx, ly, lz]),
    )
    if parse_special_pairs:
        _parse_pairs_information(gsd_snapshot, top, site_indexMap, n_rigid)
    if top.n_bonds > 0:
        _parse_bond_information(gsd_snapshot, top, site_indexMap, n_rigid)
    if top.n_angles > 0:
        _parse_angle_information(gsd_snapshot, top, site_indexMap, n_rigid)
    if top.n_dihedrals > 0:
        _parse_dihedral_information(gsd_snapshot, top, site_indexMap, n_rigid)
    if top.n_impropers > 0:
        _parse_improper_information(gsd_snapshot, top, site_indexMap, n_rigid)

    if rigid_info:
        return gsd_snapshot, base_units, rigid_info
    else:
        return gsd_snapshot, base_units


def to_hoomd_snapshot(
    top,
    base_units=None,
    shift_coords=True,
    parse_special_pairs=True,
    auto_scale=False,
):
    """Create a hoomd.Snapshot objcet (HOOMD default data format).

    The Hoomd snapshot defines the topology of HOOMD-blue simulations. This file
    can be used as a starting point for a HOOMD-blue simulation, for analysis,
    and for visualization in various tools.

    Parameters
    ----------
    top : gmso.Topology
        gmso.Topology object
    filename : str
        Path of the output file.
    base_units : dict, optinoal, default=None
        The dictionary of base units to be converted to. Entries restricted to
        "energy", "length", and "mass". There is also option to used predefined
        unit systems ("MD" or "AKMA" provided as string). If None is provided,
        this method will perform no units conversion.
    shift_coords : bool, optional, default=True
        Shift coordinates from (0, L) to (-L/2, L/2) if necessary.
    parse_special_pairs : bool, optional, default=True
        Writes out special pair information necessary to correctly use the OPLS
        fudged 1,4 interactions in HOOMD.
    auto_scale : bool or dict, optional, default=False
        Automatically scaling relevant length, energy and mass units.
        Referenced mass unit is obtained from sites' masses.
        Referenced energy and distance are refered from sites' atom types (when applicable).
        If the referenced scaling values cannot be determined (e.g., when the topology is not typed),
        all reference scaling values is set to 1.
        A dictionary specifying the referenced scaling values may also be provided for this argument.

    Return
    ------
    hoomd_snapshot : hoomd.Snapshot
        Converted hoomd Snapshot. Always retruned.
    base_units : dict
        Base units dictionary utilized during the conversion. Always returned.
    rigid_info : hoomd.md.constrain.Rigid
        A hoomd constraint object storing constituent particle information
        needed to run rigid body simulations in HOOMD-blue.
        This is only returned if sites with `site.molecule.isrigid = True`
        are found in the topology.

    Notes
    -----
    Force field parameters are not written to the GSD file and must be
    generated using `gmso.external.convert_hoomd.to_hoomd_forcefield()`.

    If you are using mBuild and GMSO to initialize rigid body simulations
    with HOOMD-blue set the `site.molecule.isrigid` property to `True`.
    If your topology contains a mix of rigid and flexible molecules,
    the rigid molecules must come first in the hierarchy of the GMSO Toplogy
    and therefore also in the mBuild Compound hierarchy.
    For more information about running rigid body simulations in HOOMD-blue,
    see https://hoomd-blue.readthedocs.io/en/latest/tutorial/06-Modelling-Rigid-Bodies/00-index.html

    For more information on the units from `base_units`, see
    https://hoomd-blue.readthedocs.io/en/latest/units.html#base-units

    Rigid Body Example
    ------------------
    In this example, a system with rigid benzene and flexible ethane
    is initialized. Since benzene is rigid, it must be passed first
    into `mb.fill_box` so that benzene molecules are at the top of
    `box` hierarchy.
        ::
            ethane = mb.lib.molecules.Ethane()
            ethane.name = "ethane"
            benzene = mb.load("c1ccccc1", smiles=True)
            benzene.name = "benzene"
            box = mb.fill_box([benzene, ethane], n_compounds=[1, 1], box=[2, 2, 2])
            top = from_mbuild(box)
            for site in top.sites:
                if site.molecule.name == "benzene":
                    site.molecule.isrigid = True

            snapshot, refs, rigid = to_hoomd_snapshot(top)
    """
    if int(hoomd_version[0]) < 4:
        raise EngineIncompatibilityError(
            "GMSO is only compatible with HOOMD-blue >= 4.0"
        )
    base_units = _validate_base_units(base_units, top, auto_scale)
    hoomd_snapshot = hoomd.Snapshot()

    # Write box information
    (lx, ly, lz, xy, xz, yz) = _prepare_box_information(top)

    lx = lx.to_value(base_units["length"])
    ly = ly.to_value(base_units["length"])
    lz = lz.to_value(base_units["length"])

    hoomd_snapshot.configuration.box = hoomd.Box(
        Lx=lx, Ly=ly, Lz=lz, xy=xy, xz=xz, yz=yz
    )

    logger.info(
        "Only writing particle, bond, angle, proper and improper dihedral information."
        "Special pairs are not currently written to GSD files",
    )

    # Create IndexMap
    site_indexMap = {id(site): i for i, site in enumerate(top.sites)}

    n_rigid, rigid_info = _parse_particle_information(
        hoomd_snapshot,
        top,
        base_units,
        shift_coords,
        u.unyt_array([lx, ly, lz]),
    )
    if parse_special_pairs:
        _parse_pairs_information(hoomd_snapshot, top, site_indexMap, n_rigid)
    if top.n_bonds > 0:
        _parse_bond_information(hoomd_snapshot, top, site_indexMap, n_rigid)
    if top.n_angles > 0:
        _parse_angle_information(hoomd_snapshot, top, site_indexMap, n_rigid)
    if top.n_dihedrals > 0:
        _parse_dihedral_information(hoomd_snapshot, top, site_indexMap, n_rigid)
    if top.n_impropers > 0:
        _parse_improper_information(hoomd_snapshot, top, site_indexMap, n_rigid)

    hoomd_snapshot.wrap()
    if rigid_info:
        return hoomd_snapshot, base_units, rigid_info
    else:
        return hoomd_snapshot, base_units


def _parse_particle_information(
    snapshot,
    top,
    base_units,
    shift_coords,
    box_lengths,
):
    """Parse site information from topology.

    Parameters
    ----------
    snapshot : gsd.hoomd.Frame or hoomd.Snapshot
        The target Snapshot object.
    top : gmso.Topology
        Topology object holding system information.
    base_units : dict
        The dictionary holding base units (mass, length, and energy)
    shift_coords : bool
        If True, shift coordinates from (0, L) to (-L/2, L/2) if neccessary.
    box_lengths : list() of length 3
        Lengths of box in x, y, z
    """
    n_sites = top.n_sites
    length_unit = base_units["length"]
    mass_unit = base_units["mass"]
    energy_unit = base_units["energy"]

    # --- Single pass over top.sites to collect all per-site data ---
    xyz_raw = np.empty((n_sites, 3), dtype=np.float64)
    masses_raw = np.empty(n_sites, dtype=np.float64)
    charges_list = np.empty(n_sites, dtype=np.float64)
    types = [None] * n_sites
    has_rigid = False

    zero_charge = (0 * u.elementary_charge).to_value(u.elementary_charge)
    default_mass = 1.0  # dimensionless, matching `1 * base_units["mass"]` via to_value

    for i, site in enumerate(top.sites):
        xyz_raw[i] = site.position.to_value(length_unit)
        types[i] = site.name if site.atom_type is None else site.atom_type.name
        masses_raw[i] = site.mass.to_value(mass_unit) if site.mass else default_mass
        charges_list[i] = (
            site.charge.to_value(u.elementary_charge) if site.charge else zero_charge
        )
        if not has_rigid and site.molecule.isrigid:
            has_rigid = True

    # Dimensionless — to_value() already converted to base unit magnitudes
    xyz = u.unyt_array(xyz_raw)
    masses = u.unyt_array(masses_raw)
    charges = u.unyt_array(charges_list, u.elementary_charge)

    if shift_coords:
        logger.info("Shifting coordinates to [-L/2, L/2]")
        xyz = coord_shift(xyz, box_lengths)

    # --- Build type IDs with dict lookup instead of list.index() ---
    unique_types = sorted(set(types))
    type_to_id = {t: idx for idx, t in enumerate(unique_types)}
    typeids = np.array([type_to_id[t] for t in types], dtype=np.int32)

    # Default orientations and moments of inertia
    moment_of_inertias = np.tile([1.0, 0.0, 0.0], (n_sites, 1))
    orientations = np.tile([1.0, 0.0, 0.0, 0.0], (n_sites, 1))

    # --- Rigid body handling ---
    if has_rigid:
        rigid_ids = np.empty(n_sites, dtype=np.int64)
        rigid_body_sets = {}

        for i, site in enumerate(top.sites):
            if site.molecule.isrigid:
                mol_num = site.molecule.number
                mol_name = site.molecule.name
                rigid_ids[i] = mol_num
                if mol_name not in rigid_body_sets:
                    rigid_body_sets[mol_name] = set()
                rigid_body_sets[mol_name].add(mol_num)
            else:
                rigid_ids[i] = -1

        # Validate hierarchy: all rigid sites must precede non-rigid sites
        non_rigid_mask = rigid_ids == -1
        if non_rigid_mask.any():
            first_neg_index = int(np.argmax(non_rigid_mask))
            if not non_rigid_mask[first_neg_index:].all():
                raise RuntimeError(
                    "When using a combination of rigid and non-rigid molecules, "
                    "all of the rigid molecules must come first in the mBuild/GMSO hierarchy."
                )

        n_rigid_types = len(rigid_body_sets)
        n_rigid = sum(len(mols) for mols in rigid_body_sets.values())

        # Dimensionless placeholders to match xyz/masses
        rigid_charges = np.zeros(n_rigid) * charges.units
        rigid_masses = u.unyt_array(np.zeros(n_rigid))
        rigid_xyz = u.unyt_array(np.zeros((n_rigid, 3)))
        rigid_moits = u.unyt_array(np.zeros((n_rigid, 3)))
        rigid_orientations = np.tile([1.0, 0.0, 0.0, 0.0], (n_rigid, 1))

        # Build rigid type IDs
        rigid_type_ids = np.empty(n_rigid, dtype=np.int32)
        offset = 0
        for i, rigid_type in enumerate(rigid_body_sets):
            count = len(rigid_body_sets[rigid_type])
            rigid_type_ids[offset : offset + count] = i
            offset += count

        # Precompute group indices for each rigid body ID (avoids repeated np.where)
        group_indices_map = {}
        for i in range(n_sites):
            rid = rigid_ids[i]
            if rid != -1:
                if rid not in group_indices_map:
                    group_indices_map[rid] = []
                group_indices_map[rid].append(i)
        for rid in group_indices_map:
            group_indices_map[rid] = np.array(group_indices_map[rid], dtype=np.intp)

        rigid_constraint = hoomd.md.constrain.Rigid()
        mol_count = 0

        for rigid_mol in rigid_body_sets:
            sorted_ids = sorted(rigid_body_sets[rigid_mol])
            for idx, _id in enumerate(sorted_ids):
                group_indices = group_indices_map[_id]
                group_positions = xyz[group_indices]
                group_masses = masses[group_indices]

                total_mass = group_masses.sum()
                com_xyz = (group_positions.T * group_masses).sum(axis=1) / total_mass

                pos_idx = idx + mol_count
                rigid_masses[pos_idx] = total_mass
                rigid_xyz[pos_idx] = com_xyz
                rigid_moits[pos_idx] = moment_of_inertia(
                    group_positions, group_masses, com_xyz
                )

                if idx == 0:
                    const_ids = typeids[group_indices]
                    const_types = [unique_types[i] for i in const_ids]
                    group_orientations = [[1, 0, 0, 0] for _ in group_indices]
                    rigid_constraint.body[rigid_mol] = {
                        "constituent_types": const_types,
                        "positions": group_positions,
                        "orientations": group_orientations,
                    }
            mol_count += len(rigid_body_sets[rigid_mol])

        # Prepend rigid body data
        unique_types = list(rigid_body_sets.keys()) + unique_types
        typeids = np.concatenate((rigid_type_ids, typeids + n_rigid_types))
        masses = np.concatenate((rigid_masses, masses))
        xyz = np.concatenate((rigid_xyz, xyz))
        charges = np.concatenate((rigid_charges, charges))
        orientations = np.concatenate((rigid_orientations, orientations), axis=0)
        moment_of_inertias = np.concatenate((rigid_moits, moment_of_inertias), axis=0)
        rigid_id_tags = np.concatenate((np.arange(n_rigid), rigid_ids))
    else:
        n_rigid = 0
        rigid_constraint = None

    # Compute charge factor
    e0 = u.physical_constants.eps_0.in_units(
        u.elementary_charge**2 / (energy_unit * length_unit)
    )
    charge_factor = (4.0 * np.pi * e0 * length_unit * energy_unit) ** 0.5

    # --- Write snapshot (deduplicated) ---
    is_hoomd_snapshot = isinstance(snapshot, hoomd.Snapshot)
    snapshot.particles.N = n_sites + n_rigid
    snapshot.particles.types = unique_types

    if is_hoomd_snapshot:
        snapshot.particles.position[0:] = xyz
        snapshot.particles.typeid[0:] = typeids
        snapshot.particles.mass[0:] = masses
        snapshot.particles.charge[0:] = charges / charge_factor
        snapshot.particles.orientation[0:] = orientations
        if n_rigid:
            snapshot.particles.body[0:] = rigid_id_tags
            snapshot.particles.moment_inertia[0:] = moment_of_inertias
    else:  # gsd.hoomd.Frame
        snapshot.particles.position = xyz
        snapshot.particles.typeid = typeids
        snapshot.particles.mass = masses
        snapshot.particles.charge = charges / charge_factor
        snapshot.particles.orientation = orientations
        if n_rigid:
            snapshot.particles.body = rigid_id_tags
            snapshot.particles.moment_inertia = moment_of_inertias

    return n_rigid, rigid_constraint


def _parse_pairs_information(snapshot, top, site_indexMap, n_rigid=0):
    """Parse sacled pair types.

    Parameters
    ----------
    snapshot : gsd.hoomd.Frame or hoomd.Snapshot
        The target Snapshot object.
    top : gmso.Topology
        Topology object holding system information
    n_rigid : int
        The number of rigid bodies found in `_parse_particle_information()`
        Used to adjust pair group indices.

    """
    pair_types = list()
    pair_typeids = list()
    pairs = list()

    scaled_pairs = list()
    pairs_dict = generate_pairs_lists(top, refer_from_scaling_factor=True)
    for pair_type in pairs_dict:
        scaled_pairs.extend(pairs_dict[pair_type])

    for pair in scaled_pairs:
        if pair[0].atom_type and pair[1].atom_type:
            pair.sort(key=lambda site: site.atom_type.name)
            pair_type = "-".join([pair[0].atom_type.name, pair[1].atom_type.name])
        else:
            pair.sort(key=lambda site: site.name)
            pair_type = "-".join([pair[0].name, pair[1].name])
        if pair_type not in pair_types:
            pair_types.append(pair_type)
        pair_typeids.append(pair_types.index(pair_type))
        pairs.append(
            (site_indexMap[id(pair[0])] + n_rigid, site_indexMap[id(pair[1])] + n_rigid)
        )

    if isinstance(snapshot, hoomd.Snapshot):
        snapshot.pairs.N = len(pairs)
        snapshot.pairs.group[:] = np.reshape(pairs, (-1, 2))
        snapshot.pairs.types = pair_types
        snapshot.pairs.typeid[:] = pair_typeids
    elif isinstance(snapshot, gsd.hoomd.Frame):
        snapshot.pairs.N = len(pairs)
        snapshot.pairs.group = np.reshape(pairs, (-1, 2))
        snapshot.pairs.types = pair_types
        snapshot.pairs.typeid = pair_typeids


def _parse_bond_information(snapshot, top, site_indexMap, n_rigid=0):
    """Parse bonds information from topology.

    Parameters
    ----------
    snapshot : gsd.hoomd.Frame or hoomd.Snapshot
        The target Snapshot object.
    top : gmso.Topology
        Topology object holding system information
    n_rigid : int
        The number of rigid bodies found in `_parse_particle_information()`
        Used to adjust bond group indices.

    """
    snapshot.bonds.N = top.n_bonds
    logger.info(f"{top.n_bonds} bonds detected")
    bond_groups = []
    bond_typeids = []
    bond_types = []

    for bond in top.bonds:
        if all([site.atom_type for site in bond.connection_members]):
            connection_members = sort_connection_members(bond, "atomclass")
            bond_type = "-".join(
                [site.atom_type.atomclass for site in connection_members]
            )
        else:
            connection_members = sort_connection_members(bond, "name")
            bond_type = "-".join([site.name for site in connection_members])

        bond_types.append(bond_type)
        bond_groups.append(
            sorted(
                tuple(site_indexMap[id(site)] + n_rigid for site in connection_members)
            )
        )

    unique_bond_types = list(set(bond_types))
    type_to_id = {btype: i for i, btype in enumerate(unique_bond_types)}
    bond_typeids = [type_to_id[btype] for btype in bond_types]  # O(n)

    if isinstance(snapshot, hoomd.Snapshot):
        snapshot.bonds.types = unique_bond_types
        snapshot.bonds.typeid[:] = bond_typeids
        snapshot.bonds.group[:] = bond_groups
    elif isinstance(snapshot, gsd.hoomd.Frame):
        snapshot.bonds.types = unique_bond_types
        snapshot.bonds.typeid = bond_typeids
        snapshot.bonds.group = bond_groups

    logger.info(f"{len(unique_bond_types)} unique bond types detected")


def _parse_angle_information(snapshot, top, site_indexMap, n_rigid=0):
    """Parse angles information from topology.

    Parameters
    ----------
    snapshot : gsd.hoomd.Frame or hoomd.Snapshot
        The target Snapshot object.
    top : gmso.Topology
        Topology object holding system information
    n_rigid : int
        The number of rigid bodies found in `_parse_particle_information()`
        Used to adjust angle group indices.

    """
    snapshot.angles.N = top.n_angles
    unique_angle_types = set()
    angle_typeids = []
    angle_groups = []
    angle_types = []

    for angle in top.angles:
        if all([site.atom_type for site in angle.connection_members]):
            connection_members = sort_connection_members(angle, "atomclass")
            angle_type = "-".join(
                [site.atom_type.atomclass for site in connection_members]
            )
        else:
            connection_members = sort_connection_members(angle, "name")
            angle_type = "-".join([site.name for site in connection_members])

        angle_types.append(angle_type)
        angle_groups.append(
            sorted(
                tuple(site_indexMap[id(site)] + n_rigid for site in connection_members)
            )
        )

    unique_angle_types = list(set(angle_types))
    type_to_id = {atype: i for i, atype in enumerate(unique_angle_types)}
    angle_typeids = [type_to_id[atype] for atype in angle_types]

    if isinstance(snapshot, hoomd.Snapshot):
        snapshot.angles.types = unique_angle_types
        snapshot.angles.typeid[:] = angle_typeids
        snapshot.angles.group[:] = np.reshape(angle_groups, (-1, 3))
    elif isinstance(snapshot, gsd.hoomd.Frame):
        snapshot.angles.types = unique_angle_types
        snapshot.angles.typeid = angle_typeids
        snapshot.angles.group = np.reshape(angle_groups, (-1, 3))

    logger.info(f"{top.n_angles} angles detected")
    logger.info(f"{len(unique_angle_types)} unique angle types detected")


def _parse_dihedral_information(snapshot, top, site_indexMap, n_rigid=0):
    """Parse dihedral information from topology.

    Parameters
    ----------
    snapshot : gsd.hoomd.Frame or hoomd.Snapshot
        The target Snapshot object.
    top : gmso.Topology
        Topology object holding system information
    n_rigid : int
        The number of rigid bodies found in `_parse_particle_information()`
        Used to adjust dihedral group indices.

    """
    snapshot.dihedrals.N = top.n_dihedrals
    dihedral_groups = []
    dihedral_types = []

    for dihedral in top.dihedrals:
        if all([site.atom_type for site in dihedral.connection_members]):
            connection_members = sort_connection_members(dihedral, "atomclass")
            dihedral_type = "-".join(
                [site.atom_type.atomclass for site in connection_members]
            )
        else:
            connection_members = sort_connection_members(dihedral, "name")
            dihedral_type = "-".join([site.name for site in connection_members])

        dihedral_types.append(dihedral_type)
        dihedral_groups.append(
            sorted(
                tuple(site_indexMap[id(site)] + n_rigid for site in connection_members)
            )
        )

    unique_dihedral_types = list(set(dihedral_types))
    type_to_id = {dtype: i for i, dtype in enumerate(unique_dihedral_types)}
    dihedral_typeids = [type_to_id[dtype] for dtype in dihedral_types]

    if isinstance(snapshot, hoomd.Snapshot):
        snapshot.dihedrals.types = unique_dihedral_types
        snapshot.dihedrals.typeid[:] = dihedral_typeids
        snapshot.dihedrals.group[:] = np.reshape(dihedral_groups, (-1, 4))
    elif isinstance(snapshot, gsd.hoomd.Frame):
        snapshot.dihedrals.types = unique_dihedral_types
        snapshot.dihedrals.typeid = dihedral_typeids
        snapshot.dihedrals.group = np.reshape(dihedral_groups, (-1, 4))

    logger.info(f"{top.n_dihedrals} dihedrals detected")
    logger.info(f"{len(unique_dihedral_types)} unique dihedral types detected")


def _parse_improper_information(snapshot, top, site_indexMap, n_rigid=0):
    """Parse impropers information from topology.

    Parameters
    ----------
    snapshot : gsd.hoomd.Snaphot or hoomd.Snapshot
        The target Snapshot object.
    top : gmso.Topology
        Topology object holding system information
    n_rigid : int
        The number of rigid bodies found in `_parse_particle_information()`

    """
    snapshot.impropers.N = top.n_impropers
    improper_groups = []
    improper_types = []

    for improper in top.impropers:
        if all([site.atom_type for site in improper.connection_members]):
            connection_members = sort_connection_members(improper, "atomclass")
            improper_type = "-".join(
                [site.atom_type.atomclass for site in connection_members]
            )
        else:
            connection_members = sort_connection_members(improper, "name")
            improper_type = "-".join([site.name for site in connection_members])

        improper_types.append(improper_type)
        improper_groups.append(
            sorted(
                tuple(site_indexMap[id(site)] + n_rigid for site in connection_members)
            )
        )

    unique_improper_types = list(set(improper_types))
    type_to_id = {itype: i for i, itype in enumerate(unique_improper_types)}
    improper_typeids = [type_to_id[itype] for itype in improper_types]

    if isinstance(snapshot, hoomd.Snapshot):
        snapshot.impropers.types = unique_improper_types
        snapshot.impropers.typeid[0:] = improper_typeids
        snapshot.impropers.group[0:] = np.reshape(improper_groups, (-1, 4))
    elif isinstance(snapshot, gsd.hoomd.Frame):
        snapshot.impropers.types = unique_improper_types
        snapshot.impropers.typeid = improper_typeids
        snapshot.impropers.group = np.reshape(improper_groups, (-1, 4))

    logger.info(f"{top.n_impropers} impropers detected")
    logger.info(f"{len(unique_improper_types)} unique dihedral types detected")


def _prepare_box_information(top):
    """Prepare the box information for writing to gsd."""
    if allclose_units(
        top.box.angles, np.array([90, 90, 90]) * u.degree, rtol=1e-5, atol=1e-8
    ):
        logger.info("Orthorhombic box detected")
        lx, ly, lz = top.box.lengths
        xy, xz, yz = 0.0, 0.0, 0.0
    else:
        logger.info("Non-orthorhombic box detected")
        u_vectors = top.box.get_unit_vectors()
        lx, ly, lz = top.box.lengths
        xy = u_vectors[1][0]
        xz = u_vectors[2][0]
        yz = u_vectors[2][1]
    return lx, ly, lz, xy, xz, yz


def to_hoomd_forcefield(
    top,
    r_cut,
    nlist=None,
    pppm_kwargs={"resolution": (8, 8, 8), "order": 4},
    base_units=None,
    auto_scale=False,
):
    """Convert the potential portion of a typed GMSO to hoomd forces.

    Parameters
    ----------
    top : gmso.Topology
        The typed topology to be converted
    r_cut : float
        r_cut for the nonbonded forces.
    nlist : hoomd.md.nlist.NeighborList or tuple or list or None, optional, default=None
        Neighborlist object to use for nonbonded forcefield. Can also be a list or tuple of neighborlists
        where the first neighborlist is used for nonbonded pair forces, and the second with coulombic forces.
        If None, the default value used will be a hoomd.md.nlist.Cell(exclusions=exclusions, buffer=0.4).
    pppm_kwargs : dict
        Keyword arguments to pass to hoomd.md.long_range.make_pppm_coulomb_forces().
    base_units : dict or str, optional, default=None
        The dictionary of base units to be converted to. Entries restricted to
        "energy", "length", and "mass". There is also option to used predefined
        unit systems ("MD" or "AKMA" provided as string). If None is provided,
        this method will perform no units conversion.
    auto_scale : bool or dict, optional, default=False
        Automatically scaling relevant length, energy and mass units.
        Referenced mass unit is obtained from sites' masses.
        Referenced energy and distance are refered from sites' atom types (when applicable).
        If the referenced scaling values cannot be determined (e.g., when the topology is not typed),
        all reference scaling values is set to 1.
        A dictionary specifying the referenced scaling values may also be provided for this argument.

    Returns
    -------
    forces : dict
        HOOMD forces converted from all available PotentialTypes of the provided
        GMSO Topology. Converted are grouped by their category (as key of the
        dictionary), namely, "nonbonded", "bonds", rangles", "dihedrals", and "impropers".
    base_units : dict
        Base units dictionary utilized during the conversion.

    """
    if int(hoomd_version[0]) < 4:
        raise EngineIncompatibilityError(
            "GMSO is only compatible with HOOMD-blue >= 4.0"
        )
    potential_types = _validate_compatibility(top)
    base_units = _validate_base_units(base_units, top, auto_scale, potential_types)

    # Reference json dict of all the potential in the PotentialTemplate
    potential_refs = dict()
    for json_file in PotentialTemplateLibrary().json_refs:
        with open(json_file) as f:
            cont = json.load(f)
        potential_refs[cont["name"]] = cont

    # convert nonbonded potentials
    forces = {
        "nonbonded": _parse_nonbonded_forces(
            top,
            r_cut,
            nlist,
            potential_types,
            potential_refs,
            pppm_kwargs,
            base_units,
        ),
        "bonds": _parse_bond_forces(
            top,
            potential_types,
            potential_refs,
            base_units,
        ),
        "angles": _parse_angle_forces(
            top,
            potential_types,
            potential_refs,
            base_units,
        ),
        "dihedrals": _parse_dihedral_forces(
            top,
            potential_types,
            potential_refs,
            base_units,
        ),
        "impropers": _parse_improper_forces(
            top,
            potential_types,
            potential_refs,
            base_units,
        ),
    }

    return forces, base_units


def _validate_compatibility(top):
    """Check and sort all the potential objects in the topology."""
    from gmso.utils.compatibility import check_compatibility

    templates = PotentialTemplateLibrary()
    lennard_jones_potential = templates["LennardJonesPotential"]
    harmonic_bond_potential = templates["HarmonicBondPotential"]
    fene_bond_potential = templates["HOOMDFENEWCABondPotential"]
    harmonic_angle_potential = templates["HarmonicAnglePotential"]
    periodic_torsion_potential = templates["PeriodicTorsionPotential"]
    hoomd_periodic_torsion_potential = templates["HOOMDPeriodicDihedralPotential"]
    harmonic_torsion_potential = templates["HarmonicTorsionPotential"]
    opls_torsion_potential = templates["OPLSTorsionPotential"]
    rb_torsion_potential = templates["RyckaertBellemansTorsionPotential"]
    accepted_potentials = (
        lennard_jones_potential,
        harmonic_bond_potential,
        fene_bond_potential,
        harmonic_angle_potential,
        periodic_torsion_potential,
        hoomd_periodic_torsion_potential,
        harmonic_torsion_potential,
        opls_torsion_potential,
        rb_torsion_potential,
    )
    potential_types = check_compatibility(top, accepted_potentials)
    return potential_types


def _parse_nonbonded_forces(
    top,
    r_cut,
    nlist,
    potential_types,
    potential_refs,
    pppm_kwargs,
    base_units,
):
    """Parse nonbonded forces from topology.

    Parameters
    ----------
    top : gmso.Topology
        Topology object holding system information.
    r_cut : float
        Cut-off radius in simulation units
    nlist : hoomd.md.nlist.NeighborList or tuple or list
        neighborlist object to use for nonbonded forcefield. Can also be a list or tuple of neighborlists
        where the first neighborlist is used for nonbonded pair forces, and the second with coulombic forces.
    potential_types : dict
        Output from _validate_compatibility().
    potential_refs : dict
        Reference json potential from gmso.lib.potential_templates.
    pppm_kwargs : dict
        Keyword arguments to pass to hoomd.md.long_range.make_pppm_coulomb_forces().
    base_units : dict
        The dictionary holding base units (mass, length, and energy)
    """
    unique_atypes = top.atom_types(filter_by=PotentialFilters.UNIQUE_NAME_CLASS)

    # Grouping atomtype by group name
    groups = dict()
    for atype in unique_atypes:
        group = potential_types[atype]
        if group not in groups:
            groups[group] = [atype]
        else:
            groups[group].append(atype)

    # Perform units conversion based on the provided base_units
    for group in groups:
        expected_units_dim = potential_refs[group]["expected_parameters_dimensions"]
        groups[group] = convert_params_units(
            groups[group],
            expected_units_dim,
            base_units,
        )

    atype_parsers = {
        "LennardJonesPotential": _parse_lj,
        "BuckinghamPotential": _parse_buckingham,
        "MiePotential": _parse_mie,
    }

    # Use Topology scaling factor to determine exclusion
    # TODO: Use molecule scaling factor
    nb_scalings, coulombic_scalings = top.scaling_factors
    if nlist is None:
        nlist_nb, nlist_coul = get_cell_nlist(top)
    elif isinstance(nlist, hoomd.md.nlist.NeighborList):
        nlist_nb = nlist
        nlist_coul = nlist
    elif isinstance(nlist, (tuple, list)):
        nlist_nb = nlist[0]
        nlist_coul = nlist[1]
    else:
        raise ValueError(
            "Incorrect values supplied for nlist. Should be of type hoomd.md.nlist"
        )

    nbonded_forces = list()
    nbonded_forces.extend(
        _parse_coulombic(
            top=top,
            nlist=nlist_coul,
            scaling_factors=coulombic_scalings,
            resolution=pppm_kwargs["resolution"],
            order=pppm_kwargs["order"],
            r_cut=r_cut,
        )
    )
    for group in groups:
        nbonded_forces.extend(
            atype_parsers[group](
                top=top,
                atypes=groups[group],
                combining_rule=top.combining_rule,
                r_cut=r_cut,
                nlist=nlist_nb,
                scaling_factors=nb_scalings,
            )
        )

    return nbonded_forces


def _parse_coulombic(
    top,
    nlist,
    scaling_factors,
    resolution,
    order,
    r_cut,
):
    """Parse coulombic forces."""
    charge_groups = any(
        [site.charge.to_value(u.elementary_charge) for site in top.sites]
    )
    if not charge_groups:
        logger.info("No charged group detected, skipping electrostatics.")
        return []
    else:
        coulombic = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(
            nlist=nlist, resolution=resolution, order=order, r_cut=r_cut
        )
    if not np.any(scaling_factors):
        return [*coulombic]  # return early
    # Handle 1-2, 1-3, and 1-4 scaling
    # TODO: Figure out a more general way to do this and handle molecule scaling factors
    special_coulombic = hoomd.md.special_pair.Coulomb()

    # Use same method as to_hoomd_snapshot to generate pairs list
    pairs_dict = generate_pairs_lists(top)
    for i, pair_type in enumerate(pairs_dict):
        if scaling_factors[i] and pairs_dict[pair_type]:
            for pair in pairs_dict[pair_type]:
                pair_name = "-".join(
                    sorted([pair[0].atom_type.name, pair[1].atom_type.name])
                )
                special_coulombic.params[pair_name] = dict(alpha=scaling_factors[i])
                special_coulombic.r_cut[pair_name] = r_cut
    return [*coulombic, special_coulombic]


def _parse_lj(top, atypes, combining_rule, r_cut, nlist, scaling_factors):
    """Parse LJ forces and special pairs LJ forces."""
    lj = hoomd.md.pair.LJ(nlist=nlist)
    calculated_params = dict()
    for pairs in itertools.combinations_with_replacement(atypes, 2):
        pairs = list(pairs)
        pairs.sort(key=lambda atype: atype.name)
        type_name = (pairs[0].name, pairs[1].name)
        comb_epsilon = np.sqrt(
            pairs[0].parameters["epsilon"].value * pairs[1].parameters["epsilon"].value
        )
        if top.combining_rule == "lorentz":
            comb_sigma = np.mean(
                [pairs[0].parameters["sigma"], pairs[1].parameters["sigma"]]
            )
        elif top.combining_rule == "geometric":
            comb_sigma = np.sqrt(
                pairs[0].parameters["sigma"].value * pairs[1].parameters["sigma"].value
            )
        else:
            raise ValueError(f"Invalid combining rule provided ({combining_rule})")

        calculated_params[type_name] = {
            "sigma": comb_sigma,
            "epsilon": comb_epsilon,
        }
        lj.params[type_name] = calculated_params[type_name]
        lj.r_cut[(type_name)] = r_cut

    # Handle 1-2, 1-3, and 1-4 scaling
    # TODO: Figure out a more general way to do this
    # and handle molecule scaling factors
    if not np.any(scaling_factors):
        return [lj]
    special_lj = hoomd.md.special_pair.LJ()

    pairs_dict = generate_pairs_lists(top)
    for i, pair_type in enumerate(pairs_dict):
        if scaling_factors[i] and pairs_dict[pair_type]:
            for pair in pairs_dict[pair_type]:
                if pair[0].atom_type in atypes and pair[1].atom_type in atypes:
                    adjscale = scaling_factors[i]
                    pair.sort(key=lambda site: site.atom_type.name)
                    pair_name = (
                        pair[0].atom_type.name,
                        pair[1].atom_type.name,
                    )
                    scaled_epsilon = adjscale * calculated_params[pair_name]["epsilon"]
                    sigma = calculated_params[pair_name]["sigma"]
                    special_lj.params["-".join(pair_name)] = {
                        "sigma": sigma,
                        "epsilon": scaled_epsilon,
                    }
                    special_lj.r_cut["-".join(pair_name)] = r_cut

    return [lj, special_lj]


# TODO: adding supports for the following nonbonded potentials
def _parse_buckingham(
    top,
    atypes,
    combining_rule,
    r_cut,
    nlist,
    scaling_factors,
):
    return None


def _parse_lj0804(
    top,
    atypes,
    combining_rule,
    r_cut,
    nlist,
    scaling_factors,
):
    return None


def _parse_lj1208(
    top,
    atypes,
    combining_rule,
    r_cut,
    nlist,
    scaling_factors,
):
    return None


def _parse_mie(
    top,
    atypes,
    combining_rule,
    r_cut,
    nlist,
    scaling_factors,
):
    return None


def _parse_bond_forces(
    top,
    potential_types,
    potential_refs,
    base_units,
):
    """Parse bond forces from topology.

    Parameters
    ----------
    top : gmso.Topology
        Topology object holding system information
    potential_types : dict
        Output from _validate_compatibility().
    potential_refs : dict
        Reference json potential from gmso.lib.potential_templates.
    base_units : dict
        The dictionary holding base units (mass, length, and energy)
    """
    unique_btypes = top.bond_types(filter_by=PotentialFilters.UNIQUE_NAME_CLASS)
    groups = dict()
    for btype in unique_btypes:
        group = potential_types[btype]
        if group not in groups:
            groups[group] = [btype]
        else:
            groups[group].append(btype)

    for group in groups:
        expected_units_dim = potential_refs[group]["expected_parameters_dimensions"]
        groups[group] = convert_params_units(
            groups[group],
            expected_units_dim,
            base_units,
        )

    btype_group_map = {
        "HarmonicBondPotential": {
            "container": hoomd.md.bond.Harmonic,
            "parser": _parse_harmonic_bond,
        },
        "HOOMDFENEWCABondPotential": {
            "container": hoomd.md.bond.FENEWCA,
            "parser": _parse_fene_bond,
        },
    }
    bond_forces = list()
    for group in groups:
        bond_forces.append(
            btype_group_map[group]["parser"](
                container=btype_group_map[group]["container"](),
                btypes=groups[group],
            )
        )
    return bond_forces


def _parse_harmonic_bond(
    container,
    btypes,
):
    for btype in btypes:
        members = sort_by_classes(btype)
        # If wild card in class, sort by types instead
        if "*" in members:
            members = sort_by_types(btype)
        container.params["-".join(members)] = {
            "k": btype.parameters["k"],
            "r0": btype.parameters["r_eq"],
        }
    return container


def _parse_fene_bond(
    container,
    btypes,
):
    for btype in btypes:
        members = sort_by_classes(btype)
        # If wild card in class, sort by types instead
        if "*" in members:
            members = sort_by_types(btype)
        container.params["-".join(members)] = {
            key: btype.parameters.get(key, 0.0)  # defaults to 0 if no sigma or epsilon
            for key in ["k", "r0", "epsilon", "sigma", "delta"]
        }
    return container


def _parse_angle_forces(
    top,
    potential_types,
    potential_refs,
    base_units,
):
    """Parse angle forces from topology.

    Parameters
    ----------
    top : gmso.Topology
        Topology object holding system information
    potential_types : dict
        Output from _validate_compatibility().
    potential_refs : dict
        Reference json potential from gmso.lib.potential_templates.
    base_units : dict
        The dictionary holding base units (mass, length, and energy)
    """
    unique_agtypes = top.angle_types(filter_by=PotentialFilters.UNIQUE_NAME_CLASS)
    groups = dict()
    for agtype in unique_agtypes:
        group = potential_types[agtype]
        if group not in groups:
            groups[group] = [agtype]
        else:
            groups[group].append(agtype)

    for group in groups:
        expected_units_dim = potential_refs[group]["expected_parameters_dimensions"]
        groups[group] = convert_params_units(
            groups[group],
            expected_units_dim,
            base_units,
        )

    agtype_group_map = {
        "HarmonicAnglePotential": {
            "container": hoomd.md.angle.Harmonic,
            "parser": _parse_harmonic_angle,
        },
    }
    angle_forces = list()
    for group in groups:
        angle_forces.append(
            agtype_group_map[group]["parser"](
                container=agtype_group_map[group]["container"](),
                agtypes=groups[group],
            )
        )
    return angle_forces


def _parse_harmonic_angle(
    container,
    agtypes,
):
    for agtype in agtypes:
        members = sort_by_classes(agtype)
        # If wild card in class, sort by types instead
        if "*" in members:
            members = sort_by_types(agtype)
        container.params["-".join(members)] = {
            "k": agtype.parameters["k"],
            "t0": agtype.parameters["theta_eq"],
        }
    return container


def _parse_dihedral_forces(
    top,
    potential_types,
    potential_refs,
    base_units,
):
    """Parse dihedral forces from topology.

    Parameters
    ----------
    top : gmso.Topology
        Topology object holding system information
    potential_types : dict
        Output from _validate_compatibility().
    potential_refs : dict
        Reference json potential from gmso.lib.potential_templates.
    base_units : dict
        The dictionary holding base units (mass, length, and energy)
    """
    unique_dihedrals = {}
    for dihedral in top.dihedrals:
        unique_members = tuple(
            [site.atom_type.atomclass for site in dihedral.connection_members]
        )
        unique_dihedrals[unique_members] = dihedral
    groups = dict()
    for dihedral in unique_dihedrals.values():
        group = potential_types[dihedral.dihedral_type]
        if group not in groups:
            groups[group] = [dihedral]
        else:
            groups[group].append(dihedral)

    expected_unitsDict = {}
    for group in groups:
        expected_unitsDict[group] = potential_refs[group][
            "expected_parameters_dimensions"
        ]
    dtype_group_map = {
        "OPLSTorsionPotential": {
            "container": hoomd.md.dihedral.OPLS,
            "parser": _parse_opls_dihedral,
        },
        "RyckaertBellemansTorsionPotential": {
            "container": hoomd.md.dihedral.OPLS,  # RBTorsion will converted to OPLS
            "parser": _parse_rb_dihedral,
        },
    }

    if float(f"{hoomd_version[0]}.{hoomd_version[1]}") >= 4.5:
        dtype_group_map["PeriodicTorsionPotential"] = {
            "container": hoomd.md.dihedral.Periodic,
            "parser": _parse_periodic_dihedral,
        }

    else:  # Hoomd < 4.5
        dtype_group_map["PeriodicTorsionPotential"] = {
            "container": hoomd.md.dihedral.Harmonic,
            "parser": _parse_periodic_dihedral,
        }
        dtype_group_map["HOOMDPeriodicTorsionPotential"] = {
            "container": hoomd.md.dihedral.Harmonic,
            "parser": _parse_hoomd_periodic_dihedral,
        }

    dihedral_forces = list()
    for group in groups:
        container = dtype_group_map[group]["container"]
        if isinstance(container(), hoomd.md.dihedral.OPLS):
            dihedral_forces.append(
                dtype_group_map[group]["parser"](
                    container=container(),
                    dihedrals=groups[group],
                    expected_units_dim=expected_unitsDict[group],
                    base_units=base_units,
                )
            )
        elif isinstance(container(), hoomd.md.dihedral.Periodic):
            dihedral_forces.extend(
                dtype_group_map[group]["parser"](
                    container=dtype_group_map[group]["container"](),
                    dihedrals=groups[group],
                    expected_units_dim=expected_unitsDict[group],
                    base_units=base_units,
                )
            )
    return dihedral_forces


def _parse_periodic_dihedral(container, dihedrals, expected_units_dim, base_units):
    containersList = []
    for _ in range(5):
        containersList.append(copy.deepcopy(container))
    for dihedral in dihedrals:
        dtype = dihedral.dihedral_type
        dtype = _convert_single_param_units(dtype, expected_units_dim, base_units)
        member_sites = sort_connection_members(dihedral, "atomclass")
        member_classes = [site.atom_type.atomclass for site in member_sites]
        if isinstance(dtype.parameters["k"], u.array.unyt_quantity):
            containersList[0].params["-".join(member_classes)] = {
                "k": dtype.parameters["k"].to_value() * 2,  # convert to K/2 form
                "d": 1,
                "n": dtype.parameters["n"].to_value(),
                "phi0": dtype.parameters["phi_eq"].to_value(),
            }
        elif isinstance(dtype.parameters["k"], u.array.unyt_array):
            paramsLen = len(dtype.parameters["k"])
            for nIndex in range(paramsLen):
                containersList[nIndex].params["-".join(member_classes)] = {
                    "k": dtype.parameters["k"].to_value()[nIndex] * 2,
                    "d": 1,
                    "n": dtype.parameters["n"].to_value()[nIndex],
                    "phi0": dtype.parameters["phi_eq"].to_value()[nIndex],
                }
    filled_containersList = []
    for i in range(5):  # take only periodic terms that have parameters
        if len(tuple(containersList[i].params.keys())) == 0:
            continue
        # add in extra parameters
        for key in containersList[0].params.keys():
            if key not in tuple(containersList[i].params.keys()):
                containersList[i].params[key] = {
                    "k": 0,
                    "d": 1,
                    "n": 0,
                    "phi0": 0,
                }
        filled_containersList.append(containersList[i])
    return filled_containersList


def _parse_hoomd_periodic_dihedral(
    container, dihedrals, expected_units_dim, base_units
):
    containersList = []
    for _ in range(5):
        containersList.append(copy.deepcopy(container))
    for dihedral in dihedrals:
        dtype = dihedral.dihedral_type
        dtype = _convert_single_param_units(dtype, expected_units_dim, base_units)
        member_sites = sort_connection_members(dihedral, "atomclass")
        member_classes = [site.atom_type.atomclass for site in member_sites]
        if isinstance(dtype.parameters["k"], u.array.unyt_quantity):
            containersList[0].params["-".join(member_classes)] = {
                "k": dtype.parameters["k"].to_value(),
                "d": dtype.parameters["k"].to_value(),
                "n": dtype.parameters["n"].to_value(),
                "phi0": dtype.parameters["phi_eq"].to_value(),
            }
        elif isinstance(dtype.parameters["k"], u.array.unyt_array):
            paramsLen = len(dtype.parameters["k"])
            for nIndex in range(paramsLen):
                containersList[nIndex].params["-".join(member_classes)] = {
                    "k": dtype.parameters["k"].to_value()[nIndex],
                    "d": dtype.parameters["k"].to_value()[nIndex],
                    "n": dtype.parameters["n"].to_value()[nIndex],
                    "phi0": dtype.parameters["phi_eq"].to_value()[nIndex],
                }
    filled_containersList = []
    for i in range(5):  # take only periodic terms that have parameters
        if len(tuple(containersList[i].params.keys())) == 0:
            continue
        # add in extra parameters
        for key in containersList[0].params.keys():
            if key not in tuple(containersList[i].params.keys()):
                containersList[i].params[key] = {
                    "k": 0,
                    "d": 1,
                    "n": 0,
                    "phi0": 0,
                }
        filled_containersList.append(containersList[i])
    return filled_containersList


def _parse_opls_dihedral(container, dihedrals, expected_units_dim, base_units):
    for dihedral in dihedrals:
        dtype = dihedral.dihedral_type
        dtype = _convert_single_param_units(dtype, expected_units_dim, base_units)
        # TODO: The range of ks is mismatched (GMSO go from k0 to k5)
        # May need to do a check that k0 == k5 == 0 or raise a warning
        member_sites = sort_connection_members(dihedral, "atomclass")
        member_classes = [site.atom_type.atomclass for site in member_sites]
        container.params["-".join(member_classes)] = {
            "k1": dtype.parameters["k1"],
            "k2": dtype.parameters["k2"],
            "k3": dtype.parameters["k3"],
            "k4": dtype.parameters["k4"],
        }
    return container


def _parse_rb_dihedral(container, dihedrals, expected_units_dim, base_units):
    logger.info(
        "RyckaertBellemansTorsionPotential will be converted to OPLSTorsionPotential."
    )
    for dihedral in dihedrals:
        dtype = dihedral.dihedral_type
        dtype = _convert_single_param_units(dtype, expected_units_dim, base_units)
        opls = convert_ryckaert_to_opls(dtype)
        member_sites = sort_connection_members(dihedral, "atomclass")
        member_classes = [site.atom_type.atomclass for site in member_sites]
        # TODO: The range of ks is mismatched (GMSO go from k0 to k5)
        # May need to do a check that k0 == k5 == 0 or raise a warning
        container.params["-".join(member_classes)] = {
            "k1": opls.parameters["k1"],
            "k2": opls.parameters["k2"],
            "k3": opls.parameters["k3"],
            "k4": opls.parameters["k4"],
        }
    return container


def _parse_improper_forces(
    top,
    potential_types,
    potential_refs,
    base_units,
):
    """Parse improper forces from topology.

    Parameters
    ----------
    top : gmso.Topology
        Topology object holding system information
    potential_types : dict
        Output from _validate_compatibility().
    potential_refs : dict
        Reference json potential from gmso.lib.potential_templates.
    base_units : dict
        The dictionary holding base units (mass, length, and energy)
    """
    unique_dtypes = top.improper_types(filter_by=PotentialFilters.UNIQUE_NAME_CLASS)
    groups = dict()
    for itype in unique_dtypes:
        group = potential_types[itype]
        if group not in groups:
            groups[group] = [itype]
        else:
            groups[group].append(itype)

    for group in groups:
        expected_units_dim = potential_refs[group]["expected_parameters_dimensions"]
        groups[group] = convert_params_units(
            groups[group],
            expected_units_dim,
            base_units,
        )

    if float(f"{hoomd_version[0]}.{hoomd_version[1]}") >= 4.5:
        itype_group_map = {
            "HarmonicImproperPotential": {
                "container": hoomd.md.improper.Harmonic,
                "parser": _parse_harmonic_improper,
            },
            "PeriodicTorsionPotential": {
                "container": hoomd.md.improper.Periodic,
                "parser": _parse_periodic_improper,
            },
            "HOOMDPeriodicDihedralPotential": {
                "container": hoomd.md.improper.Periodic,
                "parser": _parse_hoomd_periodic_improper,
            },
        }
    else:
        itype_group_map = {
            "HarmonicImproperPotential": {
                "container": hoomd.md.improper.Harmonic,
                "parser": _parse_harmonic_improper,
            },
        }

    improper_forces = list()
    for group in groups:
        improper_forces.append(
            itype_group_map[group]["parser"](
                container=itype_group_map[group]["container"](),
                itypes=groups[group],
            )
        )
    return improper_forces


def _parse_harmonic_improper(
    container,
    itypes,
):
    for itype in itypes:
        members = sort_by_classes(itype)
        # If wild card in class, sort by types instead
        if "*" in members:
            members = sort_by_types(itype)
        container.params["-".join(members)] = {
            "k": itype.parameters["k"],
            "chi0": itype.parameters["phi_eq"],
        }
    return container


def _parse_periodic_improper(
    container,
    itypes,
):
    for itype in itypes:
        members = sort_by_classes(itype)
        # If wild card in class, sort by types instead
        if "*" in members:
            members = sort_by_types(itype)
        container.params["-".join(members)] = {
            "k": itype.parameters["k"] * 2,  # convert to k/2
            "chi0": itype.parameters["phi_eq"],
            "n": itype.parameters["n"],
            "d": itype.parameters.get("d", 1.0),
        }
    return container


def _parse_hoomd_periodic_improper(
    container,
    itypes,
):
    for itype in itypes:
        members = sort_by_classes(itype)
        # If wild card in class, sort by types instead
        if "*" in members:
            members = sort_by_types(itype)
        container.params["-".join(members)] = {
            "k": itype.parameters["k"],
            "chi0": itype.parameters["phi0"],
            "n": itype.parameters["n"],
            "d": itype.parameters.get("d"),
        }
    return container


def _validate_base_units(base_units, top, auto_scale, potential_types=None):
    """Validate the provided base units, infer units (based on top's positions and masses) if none is provided."""
    if base_units and auto_scale:
        logger.info(
            "Both base_units and auto_scale are provided, auto_scale will take precedent."
        )

    base_units = copy.deepcopy(base_units)
    ref = {
        "energy": u.dimensions.energy,
        "length": u.dimensions.length,
        "mass": u.dimensions.mass,
    }
    unit_systems = {"MD": MD_UNITS, "AKMA": AKMA_UNITS}

    if auto_scale:
        base_units = _infer_units(top)

        # Refer masses from sites' masses
        masses = [site.mass.to(base_units["mass"]) for site in top.sites]
        if masses:
            base_units["mass"] = max(masses)
            # Refer lengths and energies from sites' atom types if possible
        unique_atypes = top.atom_types(
            filter_by=PotentialFilters.UNIQUE_NAME_CLASS,
        )
        if unique_atypes:
            if not potential_types:
                potential_types = _validate_compatibility(top)
            atype_classes = dict()
            # Separate atypes by their classes
            for atype in unique_atypes:
                if potential_types[atype] not in atype_classes:
                    atype_classes[potential_types[atype]] = [atype]
                else:
                    atype_classes[potential_types[atype]].append(atype)

        # Appending lengths and energy
        lengths, energies = list(), list()
        for atype_class in atype_classes:
            if atype_class == "LennardJonesPotential":
                for atype in unique_atypes:
                    lengths.append(atype.parameters["sigma"].to(base_units["length"]))
                    energies.append(
                        atype.parameters["epsilon"].to(base_units["energy"])
                    )
            else:
                raise NotYetImplementedWarning(
                    f"Currently cannot infer referenced lengths and energies from {atype_class}"
                )
        base_units["length"] = max(lengths)
        base_units["energy"] = max(energies)

    elif isinstance(base_units, str):
        base_units = unit_systems[base_units]
    elif isinstance(base_units, dict):
        for key in base_units:
            if key not in ["energy", "mass", "length"]:
                logger.info(
                    "Only base unit will be used during the conversion "
                    "i.e., energy, mass, and length, other units provided "
                    "will not be considered."
                )
            msg = "{key} is in wrong unit dimension"
            if isinstance(base_units[key], u.Unit):
                assert base_units[key].dimensions == ref[key], msg
                base_units[key] = 1 * base_units[key]
            elif isinstance(base_units[key], u.array.unyt_quantity):
                assert base_units[key].units.dimensions == ref[key], msg
            else:
                raise TypeError(
                    f"Base unit of {key} must be of type u.Unit or u.unyt_quantity."
                )

        missing = list()
        for base in ["energy", "mass", "length"]:
            if base not in base_units:
                missing.append(base)
        if missing:
            raise (f"base_units is not fully provided, missing {missing}")
    else:
        base_units = _infer_units(top)
        for key in base_units:
            if isinstance(base_units[key], u.Unit):
                base_units[key] = 1 * base_units[key]
        logger.warning(
            "Neither base_units or auto_scale is provided, "
            f"so default units of {base_units} are inferred."
        )
    # Add angle unit (since HOOMD will use radian across the board)
    base_units["angle"] = 1 * u.radian
    # add dimensionless handling
    base_units["dimensionless"] = 1 * u.dimensionless

    return base_units


def _infer_units(top):
    """Try to infer unit from topology."""
    mass_unit = u.unyt_array([site.mass for site in top.sites]).units
    length_unit = u.unyt_array([site.position for site in top.sites]).units

    if length_unit == u.angstrom:
        energy_unit = u.kcal / u.mol
    elif length_unit == u.nm:
        energy_unit = u.kJ / u.mol
    else:
        raise ValueError(f"Cannot infer energy unit from {length_unit}")

    return {"length": length_unit, "energy": energy_unit, "mass": mass_unit}


def _convert_single_param_units(
    potential,
    expected_units_dim,
    base_units,
):
    """Convert parameters' units in the potential to that specified in the base_units."""
    converted_params = dict()
    for parameter in potential.parameters:
        unit_dim = expected_units_dim[parameter]
        ind_units = re.sub("[^a-zA-Z]+", " ", unit_dim).split()
        for unit in ind_units:
            unit_dim = unit_dim.replace(
                unit,
                f"({str(base_units[unit].value)} * {str(base_units[unit].units)})",
            )

        converted_params[parameter] = potential.parameters[parameter].to(unit_dim)
    potential.parameters = converted_params
    return potential
