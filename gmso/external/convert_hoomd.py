"""Convert GMSO Topology to GSD snapshot."""
from __future__ import division

import itertools
import json
import re
import statistics
import warnings

import numpy as np
import unyt as u
from unyt.array import allclose_units

import gmso
from gmso.core.views import PotentialFilters
from gmso.exceptions import NotYetImplementedWarning
from gmso.formats.formats_registry import saves_as
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.utils.conversions import (
    convert_opls_to_ryckaert,
    convert_ryckaert_to_opls,
)
from gmso.utils.geometry import coord_shift
from gmso.utils.io import has_gsd, has_hoomd
from gmso.utils.sorting import (
    natural_sort,
    sort_connection_members,
    sort_member_types,
)

if has_gsd:
    import gsd.hoomd
if has_hoomd:
    import hoomd

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


def to_gsd_snapshot(
    top,
    base_units=None,
    rigid_bodies=None,
    shift_coords=True,
):
    """Create a gsd.snapshot objcet (HOOMD v3 default data format).

    The gsd snapshot is molecular structure of HOOMD-Blue. This file
    can be used as a starting point for a HOOMD-Blue simulation, for analysis,
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
    base_units = _validate_base_units(base_units, top)
    gsd_snapshot = gsd.hoomd.Snapshot()

    gsd_snapshot.configuration.step = 0
    gsd_snapshot.configuration.dimensions = 3

    # Write box information
    (lx, ly, lz, xy, xz, yz) = _prepare_box_information(top)

    lx = lx.to(base_units["length"])
    ly = ly.to(base_units["length"])
    lz = lz.to(base_units["length"])

    gsd_snapshot.configuration.box = np.array([lx, ly, lz, xy, xz, yz])

    warnings.warn(
        "Only writing particle, bond, angle, proper and improper dihedral information."
        "Special pairs are not currently written to GSD files",
        NotYetImplementedWarning,
    )

    _parse_particle_information(
        gsd_snapshot,
        top,
        base_units,
        rigid_bodies,
        shift_coords,
        u.unyt_array([lx, ly, lz]),
    )
    _parse_pairs_information(gsd_snapshot, top)
    if top.n_bonds > 0:
        _parse_bond_information(gsd_snapshot, top)
    if top.n_angles > 0:
        _parse_angle_information(gsd_snapshot, top)
    if top.n_dihedrals > 0:
        _parse_dihedral_information(gsd_snapshot, top)
    if top.n_impropers > 0:
        _parse_improper_information(gsd_snapshot, top)

    return gsd_snapshot


def _parse_pairs_information(
    snapshot,
    top,
):
    """Parse scaled pair types."""
    pair_types = list()
    pair_typeids = list()
    pairs = list()

    scaled_pairs = list()
    for pair_type in _generate_pairs_list(top):
        scaled_pairs.extend(pair_type)

    for pair in scaled_pairs:
        if pair[0].atom_type and pair[1].atom_type:
            pair.sort(key=lambda site: site.atom_type.name)
            pair_type = "-".join(
                [pair[0].atom_type.name, pair[1].atom_type.name]
            )
        else:
            pair.sort(key=lambda site: site.name)
            pair_type = "-".join([(pair[0].name, pair[1].name)])
        if pair_type not in pair_types:
            pair_types.append(pair_type)
        pair_typeids.append(pair_types.index(pair_type))
        pairs.append((top.get_index(pair[0]), top.get_index(pair[1])))

    if isinstance(snapshot, hoomd.Snapshot):
        snapshot.pairs.N = len(pairs)
        snapshot.pairs.group[:] = np.reshape(pairs, (-1, 2))
        snapshot.pairs.types = pair_types
        snapshot.pairs.typeid[:] = pair_typeids
    elif isinstance(snapshot, gsd.hoomd.Snapshot):
        snapshot.pairs.N = len(pairs)
        snapshot.pairs.group = np.reshape(pairs, (-1, 2))
        snapshot.pairs.types = pair_types
        snapshot.pairs.typeid = pair_typeids


def to_hoomd_snapshot(
    top,
    base_units=None,
    rigid_bodies=None,
    shift_coords=True,
):
    """Create a gsd.snapshot objcet (HOOMD v3 default data format).

    The gsd snapshot is molecular structure of HOOMD-Blue. This file
    can be used as a starting point for a HOOMD-Blue simulation, for analysis,
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
    base_units = _validate_base_units(base_units, top)
    hoomd_snapshot = hoomd.Snapshot()

    # Write box information
    (lx, ly, lz, xy, xz, yz) = _prepare_box_information(top)

    lx = lx.to(base_units["length"])
    ly = ly.to(base_units["length"])
    lz = lz.to(base_units["length"])

    hoomd_snapshot.configuration.box = hoomd.Box(
        Lx=lx, Ly=ly, Lz=lz, xy=xy, xz=xz, yz=yz
    )

    warnings.warn(
        "Only writing particle, bond, angle, proper and improper dihedral information."
        "Special pairs are not currently written to GSD files",
        NotYetImplementedWarning,
    )

    _parse_particle_information(
        hoomd_snapshot,
        top,
        base_units,
        rigid_bodies,
        shift_coords,
        u.unyt_array([lx, ly, lz]),
    )
    _parse_pairs_information(hoomd_snapshot, top)
    if top.n_bonds > 0:
        _parse_bond_information(hoomd_snapshot, top)
    if top.n_angles > 0:
        _parse_angle_information(hoomd_snapshot, top)
    if top.n_dihedrals > 0:
        _parse_dihedral_information(hoomd_snapshot, top)
    if top.n_impropers > 0:
        _parse_improper_information(hoomd_snapshot, top)

    hoomd_snapshot.wrap()
    return hoomd_snapshot


def _parse_particle_information(
    snapshot,
    top,
    base_units,
    rigid_bodies,
    shift_coords,
    box_lengths,
):
    """Parse site information from topology.

    Parameters
    ----------
    snapshot : gsd.hoomd.Snapshot or hoomd.Snapshot
        The target Snapshot object.
    top : gmso.Topology
        Topology object holding system information.
    base_units : dict
        The dictionary holding base units (mass, length, and energy)
    rigid_bodies : bool
        Flag to parse rigid bodies information, to be implemented
    shift_coords : bool
        If True, shift coordinates from (0, L) to (-L/2, L/2) if neccessary.
    box_lengths : list() of length 3
        Lengths of box in x, y, z
    """
    # Set up all require
    xyz = u.unyt_array([site.position for site in top.sites]).to(
        base_units["length"]
    )
    if shift_coords:
        warnings.warn("Shifting coordinates to [-L/2, L/2]")
        xyz = coord_shift(xyz, box_lengths)

    types = [
        site.name if site.atom_type is None else site.atom_type.name
        for site in top.sites
    ]
    unique_types = sorted(list(set(types)))
    typeids = np.array([unique_types.index(t) for t in types])
    masses = list()
    charges = list()
    for site in top.sites:
        masses.append(site.mass if site.mass else 1 * base_units["mass"])
        charges.append(site.charge if site.charge else 0 * u.elementary_charge)
    masses = u.unyt_array(masses).in_units(base_units["mass"]).value

    """
    Permittivity of free space = 2.39725e-4 e^2/((kcal/mol)(angstrom)),
    where e is the elementary charge
    """
    # e0 = u.physical_constants.eps_0.in_units(
    #    u.elementary_charge**2 / u.Unit("kcal*angstrom/mol")
    # )
    # charge_factor = (4.0 * np.pi * e0 * ref_distance * ref_energy) ** 0.5

    e0 = u.physical_constants.eps_0.in_units(
        u.elementary_charge**2 / (base_units["energy"] * base_units["length"])
    )
    charge_factor = (
        4.0 * np.pi * e0 * base_units["length"] * base_units["energy"]
    ) ** 0.5

    if isinstance(snapshot, hoomd.Snapshot):
        snapshot.particles.N = top.n_sites
        snapshot.particles.types = unique_types
        snapshot.particles.position[0:] = xyz.in_units(
            base_units["length"]
        ).value
        snapshot.particles.typeid[0:] = typeids
        snapshot.particles.mass[0:] = masses
        snapshot.particles.charge[0:] = charges / charge_factor
    elif isinstance(snapshot, gsd.hoomd.Snapshot):
        snapshot.particles.N = top.n_sites
        snapshot.particles.types = unique_types
        snapshot.particles.position = xyz.in_units(base_units["length"]).value
        snapshot.particles.typeid = typeids
        snapshot.particles.mass = masses
        snapshot.particles.charge = charges / charge_factor
    if rigid_bodies:
        warnings.warn(
            "Rigid bodies detected, but not yet implemented for GSD",
            NotYetImplementedWarning,
        )


def _parse_bond_information(snapshot, top):
    """Parse bonds information from topology.

    Parameters
    ----------
    snapshot : gsd.hoomd.Snapshot or hoomd.Snapshot
        The target Snapshot object.
    top : gmso.Topology
        Topology object holding system information

    """
    snapshot.bonds.N = top.n_bonds
    warnings.warn(f"{top.n_bonds} bonds detected")
    bond_groups = []
    bond_typeids = []
    bond_types = []

    for bond in top.bonds:
        if all([site.atom_type for site in bond.connection_members]):
            connection_members = sort_connection_members(bond, "atom_type")
            bond_type = "-".join(
                [site.atom_type.name for site in connection_members]
            )
        else:
            connection_members = sort_connection_members(bond, "name")
            bond_type = "-".join([site.name for site in connection_members])

        bond_types.append(bond_type)
        bond_groups.append(
            tuple(top.get_index(site) for site in connection_members)
        )

    unique_bond_types = list(set(bond_types))
    bond_typeids = [unique_bond_types.index(i) for i in bond_types]

    if isinstance(snapshot, hoomd.Snapshot):
        snapshot.bonds.types = unique_bond_types
        snapshot.bonds.typeid[0:] = bond_typeids
        snapshot.bonds.group[0:] = bond_groups
    elif isinstance(snapshot, gsd.hoomd.Snapshot):
        snapshot.bonds.types = unique_bond_types
        snapshot.bonds.typeid = bond_typeids
        snapshot.bonds.group = bond_groups

    warnings.warn(f"{len(unique_bond_types)} unique bond types detected")


def _parse_angle_information(snapshot, top):
    """Parse angles information from topology.

    Parameters
    ----------
    snapshot : gsd.hoomd.Snapshot or hoomd.Snapshot
        The target Snapshot object.
    top : gmso.Topology
        Topology object holding system information

    """
    snapshot.angles.N = top.n_angles
    unique_angle_types = set()
    angle_typeids = []
    angle_groups = []
    angle_types = []

    for angle in top.angles:
        if all([site.atom_type for site in angle.connection_members]):
            connection_members = sort_connection_members(angle, "atom_type")
            angle_type = "-".join(
                [site.atom_type.name for site in connection_members]
            )
        else:
            connection_members = sort_connection_members(angle, "name")
            angle_type = "-".join([site.name for site in connection_members])

        angle_types.append(angle_type)
        angle_groups.append(
            tuple(top.get_index(site) for site in connection_members)
        )

    unique_angle_types = list(set(angle_types))
    angle_typeids = [unique_angle_types.index(i) for i in angle_types]

    if isinstance(snapshot, hoomd.Snapshot):
        snapshot.angles.types = unique_angle_types
        snapshot.angles.typeid[0:] = angle_typeids
        snapshot.angles.group[0:] = np.reshape(angle_groups, (-1, 3))
    elif isinstance(snapshot, gsd.hoomd.Snapshot):
        snapshot.angles.types = unique_angle_types
        snapshot.angles.typeid = angle_typeids
        snapshot.angles.group = np.reshape(angle_groups, (-1, 3))

    warnings.warn(f"{top.n_angles} angles detected")
    warnings.warn(f"{len(unique_angle_types)} unique angle types detected")


def _parse_dihedral_information(snapshot, top):
    """Parse dihedral information from topology.

    Parameters
    ----------
    snapshot : gsd.hoomd.Snapshot or hoomd.Snapshot
        The target Snapshot object.
    top : gmso.Topology
        Topology object holding system information

    """
    snapshot.dihedrals.N = top.n_dihedrals
    dihedral_groups = []
    dihedral_types = []

    for dihedral in top.dihedrals:
        if all([site.atom_type for site in dihedral.connection_members]):
            connection_members = sort_connection_members(dihedral, "atom_type")
            dihedral_type = "-".join(
                [site.atom_type.name for site in connection_members]
            )
        else:
            connection_members = sort_connection_members(dihedral, "name")
            dihedral_type = "-".join([site.name for site in connection_members])

        dihedral_types.append(dihedral_type)
        dihedral_groups.append(
            tuple(top.get_index(site) for site in connection_members)
        )

    unique_dihedral_types = list(set(dihedral_types))
    dihedral_typeids = [unique_dihedral_types.index(i) for i in dihedral_types]

    if isinstance(snapshot, hoomd.Snapshot):
        snapshot.dihedrals.types = unique_dihedral_types
        snapshot.dihedrals.typeid[0:] = dihedral_typeids
        snapshot.dihedrals.group[0:] = np.reshape(dihedral_groups, (-1, 4))
    elif isinstance(snapshot, gsd.hoomd.Snapshot):
        snapshot.dihedrals.types = unique_dihedral_types
        snapshot.dihedrals.typeid = dihedral_typeids
        snapshot.dihedrals.group = np.reshape(dihedral_groups, (-1, 4))

    warnings.warn(f"{top.n_dihedrals} dihedrals detected")
    warnings.warn(
        f"{len(unique_dihedral_types)} unique dihedral types detected"
    )


def _parse_improper_information(snapshot, top):
    """Parse impropers information from topology.

    Parameters
    ----------
    snapshot : gsd.hoomd.Snaphot or hoomd.Snapshot
        The target Snapshot object.
    top : gmso.Topology
        Topology object holding system information

    """
    snapshot.impropers.N = top.n_impropers
    improper_groups = []
    improper_types = []

    for improper in top.impropers:
        if all([site.atom_type for site in improper.connection_members]):
            connection_members = sort_connection_members(improper, "atom_type")
            improper_type = "-".join(
                [site.atom_type.name for site in connection_members]
            )
        else:
            connection_members = sort_connection_members(improper, "name")
            improper_type = "-".join([site.name for site in connection_members])

        improper_types.append(improper_type)
        improper_groups.append(
            tuple(top.get_index(site) for site in connection_members)
        )

    unique_improper_types = list(set(improper_types))
    improper_typeids = [unique_improper_types.index(i) for i in improper_types]

    if isinstance(snapshot, hoomd.Snapshot):
        snapshot.impropers.types = unique_improper_types
        snapshot.impropers.typeid[0:] = improper_typeids
        snapshot.impropers.group[0:] = np.reshape(improper_groups, (-1, 4))
    elif isinstance(snapshot, gsd.hoomd.Snapshot):
        snapshot.impropers.types = unique_improper_types
        snapshot.impropers.typeid[0:] = improper_typeids
        snapshot.impropers.group[0:] = np.reshape(improper_groups, (-1, 4))

    warnings.warn(f"{top.n_impropers} impropers detected")
    warnings.warn(
        f"{len(unique_improper_types)} unique dihedral types detected"
    )


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


def to_hoomd_forcefield(
    top,
    nlist_buffer=0.4,
    r_cut=0,
    pppm_kwargs={"resolution": (8, 8, 8), "order": 4},
    base_units=None,
):
    """Convert the potential portion of a typed GMSO to hoomd forces.

    Parameters
    ----------
    top : gmso.Topology
        The typed topology to be converted
    nlist_buffer : float, optional, default=0.4
        Neighborlist buffer for simulation cell. Its unit is the same as that
        used to defined GMSO Topology Box.
    r_cut : float
        r_cut for the nonbonded forces.
    pppm_kwargs : dict
        Keyword arguments to pass to hoomd.md.long_range.make_pppm_coulomb_forces().
    base_units : dict or str, optional, default=None
        The dictionary of base units to be converted to. Entries restricted to
        "energy", "length", and "mass". There is also option to used predefined
        unit systems ("MD" or "AKMA" provided as string). If None is provided,
        this method will perform no units conversion.

    Returns
    -------
    forces : dict
        Each entry
    """
    potential_types = _validate_compatibility(top)
    base_units = _validate_base_units(base_units, top)

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
            nlist_buffer,
            r_cut,
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

    return forces


def _validate_compatibility(top):
    """Check and sort all the potential objects in the topology."""
    from gmso.utils.compatibility import check_compatibility

    templates = PotentialTemplateLibrary()
    lennard_jones_potential = templates["LennardJonesPotential"]
    harmonic_bond_potential = templates["HarmonicBondPotential"]
    harmonic_angle_potential = templates["HarmonicAnglePotential"]
    periodic_torsion_potential = templates["PeriodicTorsionPotential"]
    opls_torsion_potential = templates["OPLSTorsionPotential"]
    rb_torsion_potential = templates["RyckaertBellemansTorsionPotential"]
    accepted_potentials = (
        lennard_jones_potential,
        harmonic_bond_potential,
        harmonic_angle_potential,
        periodic_torsion_potential,
        opls_torsion_potential,
        rb_torsion_potential,
    )
    potential_types = check_compatibility(top, accepted_potentials)
    return potential_types


def _parse_nonbonded_forces(
    top,
    nlist_buffer,
    r_cut,
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
    nlist_buffer : float
        Buffer argument ot pass to hoomd.md.nlist.Cell.
    r_cut : float
        Cut-off radius in simulation units
    potential_types : dict
        Output from _validate_compatibility().
    potential_refs : dict
        Reference json potential from gmso.lib.potential_templates.
    pppm_kwargs : dict
        Keyword arguments to pass to hoomd.md.long_range.make_pppm_coulomb_forces().
    base_units : dict
        The dictionary holding base units (mass, length, and energy)
    """
    # Set up helper methods to parse different nonbonded forces.
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
            print("No charged group detected, skipping electrostatics.")
            return None
        else:
            coulombic = hoomd.md.long_range.pppm.make_pppm_coulomb_forces(
                nlist=nlist, resolution=resolution, order=order, r_cut=r_cut
            )

        # Handle 1-2, 1-3, and 1-4 scaling
        # TODO: Fiure out a more general way to do this and handle molecule scaling factors
        special_coulombic = hoomd.md.special_pair.Coulomb()

        # Use same method as to_hoomd_snapshot to generate pairs list
        for i, pairs in enumerate(_generate_pairs_list(top)):
            if scaling_factors[i] and pairs:
                for pair in pairs:
                    pair_name = "-".join(
                        [pair[0].atom_type.name, pair[1].atom_type.name]
                    )
                    special_coulombic.params[pair_name] = dict(
                        alpha=scaling_factors[i]
                    )
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
            comb_epsilon = statistics.geometric_mean(
                [pairs[0].parameters["epsilon"], pairs[1].parameters["epsilon"]]
            )
            if top.combining_rule == "lorentz":
                comb_sigma = np.mean(
                    [pairs[0].parameters["sigma"], pairs[1].parameters["sigma"]]
                )
            elif top.combining_rule == "geometric":
                comb_sigma = statistics.geometric_mean(
                    [pairs[0].parameters["sigma"], pairs[1].parameters["sigma"]]
                )
            else:
                raise ValueError(
                    f"Invalid combining rule provided ({combining_rule})"
                )

            calculated_params[type_name] = {
                "sigma": comb_sigma,
                "epsilon": comb_epsilon,
            }
            lj.params[type_name] = calculated_params[type_name]
            lj.r_cut[(type_name)] = r_cut

        # Handle 1-2, 1-3, and 1-4 scaling
        # TODO: Fiure out a more general way to do this and handle molecule scaling factors
        special_lj = hoomd.md.special_pair.LJ()

        for i, pairs in enumerate(_generate_pairs_list(top)):
            if scaling_factors[i] and pairs:
                for pair in pairs:
                    if (
                        pair[0].atom_type in atypes
                        and pair[1].atom_type in atypes
                    ):
                        adjscale = scaling_factors[i]
                        pair.sort(key=lambda site: site.atom_type.name)
                        pair_name = (
                            pair[0].atom_type.name,
                            pair[1].atom_type.name,
                        )
                        scaled_epsilon = (
                            adjscale * calculated_params[pair_name]["epsilon"]
                        )
                        sigma = calculated_params[pair_name]["sigma"]
                        special_lj.params["-".join(pair_name)] = {
                            "sigma": sigma,
                            "epsilon": scaled_epsilon,
                        }
                        special_lj.r_cut["-".join(pair_name)] = r_cut

        return [lj, special_lj]

    def _parse_buckingham(
        top, atypes, combining_rule, r_cut, nlist, scaling_factors
    ):
        return None

    def _parse_lj0804(
        top, atypes, combining_rule, r_cut, nlist, scaling_factors
    ):
        return None

    def _parse_lj1208(
        top, atypes, combining_rule, r_cut, nlist, scaling_factors
    ):
        return None

    def _parse_mie(top, atypes, combining_rule, r_cut, nlist, scaling_factors):
        return None

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
        expected_units_dim = potential_refs[group][
            "expected_parameters_dimensions"
        ]
        groups[group] = _convert_params_units(
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
    exclusions = list()
    for i in range(len(nb_scalings)):
        if i == 0:
            exclusions.append("bond")
        else:
            exclusions.append(f"1-{i+2}")
    nlist = hoomd.md.nlist.Cell(exclusions=exclusions, buffer=nlist_buffer)

    nbonded_forces = list()
    nbonded_forces.extend(
        _parse_coulombic(
            top=top,
            nlist=nlist,
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
                nlist=nlist,
                scaling_factors=nb_scalings,
            )
        )
    return nbonded_forces


def _parse_bond_forces(top, potential_types, potential_refs, base_units):
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

    def _parse_harmonic(container, btypes):
        for btype in btypes:
            # TODO: Unit conversion
            member_types = sort_member_types(btype)
            container.params["-".join(member_types)] = {
                "k": btype.parameters["k"],
                "r0": btype.parameters["r_eq"],
            }
        return container

    unique_btypes = top.bond_types(filter_by=PotentialFilters.UNIQUE_NAME_CLASS)
    groups = dict()
    for btype in unique_btypes:
        group = potential_types[btype]
        if group not in groups:
            groups[group] = [btype]
        else:
            groups[group].append(btype)

    for group in groups:
        expected_units_dim = potential_refs[group][
            "expected_parameters_dimensions"
        ]
        groups[group] = _convert_params_units(
            groups[group], expected_units_dim, base_units
        )

    btype_group_map = {
        "HarmonicBondPotential": {
            "container": hoomd.md.bond.Harmonic,
            "parser": _parse_harmonic,
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


def _parse_angle_forces(top, potential_types, potential_refs, base_units):
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

    def _parse_harmonic(container, agtypes):
        for agtype in agtypes:
            member_types = sort_member_types(agtype)
            container.params["-".join(member_types)] = {
                "k": agtype.parameters["k"],
                "t0": agtype.parameters["theta_eq"],
            }
        return container

    unique_agtypes = top.angle_types(
        filter_by=PotentialFilters.UNIQUE_NAME_CLASS
    )
    groups = dict()
    for agtype in unique_agtypes:
        group = potential_types[agtype]
        if group not in groups:
            groups[group] = [agtype]
        else:
            groups[group].append(agtype)

    base_units["angle"] = u.radian
    for group in groups:
        expected_units_dim = potential_refs[group][
            "expected_parameters_dimensions"
        ]
        groups[group] = _convert_params_units(
            groups[group], expected_units_dim, base_units
        )

    agtype_group_map = {
        "HarmonicAnglePotential": {
            "container": hoomd.md.angle.Harmonic,
            "parser": _parse_harmonic,
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


def _parse_dihedral_forces(top, potential_types, potential_refs, base_units):
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

    def _parse_periodic(container, dtypes):
        for dtype in dtypes:
            member_types = sort_member_types(dtype)
            container.params["-".join(member_types)] = {
                "k": dtype.parameters["k"],
                "d": 1,
                "n": dtype.parameters["n"],
                "phi0": dtype.parameters["phi_eq"],
            }
        return container

    def _parse_opls(container, dtypes):
        for dtype in dtypes:
            # TODO: The range of ks is mismatched (GMSO go from k0 to k5)
            # May need to do a check that k0 == k5 == 0 or raise a warning
            container.params["-".join(dtype.member_types)] = {
                "k1": dtype.parameters["k1"],
                "k2": dtype.parameters["k2"],
                "k3": dtype.parameters["k3"],
                "k4": dtype.parameters["k4"],
            }
        return container

    def _parse_rb(container, dtypes):
        warnings.warn(
            "RyckaertBellemansTorsionPotential will be converted to OPLSTorsionPotential."
        )
        for dtype in dtypes:
            opls = convert_ryckaert_to_opls(dtype)
            member_types = sort_member_types(dtype)
            # TODO: The range of ks is mismatched (GMSO go from k0 to k5)
            # May need to do a check that k0 == k5 == 0 or raise a warning
            container.params["-".join(member_types)] = {
                "k1": opls.parameters["k1"],
                "k2": opls.parameters["k2"],
                "k3": opls.parameters["k3"],
                "k4": opls.parameters["k4"],
            }
        return container

    unique_dtypes = top.dihedral_types(
        filter_by=PotentialFilters.UNIQUE_NAME_CLASS
    )
    groups = dict()
    for dtype in unique_dtypes:
        group = potential_types[dtype]
        if group not in groups:
            groups[group] = [dtype]
        else:
            groups[group].append(dtype)

    for group in groups:
        expected_units_dim = potential_refs[group][
            "expected_parameters_dimensions"
        ]
        groups[group] = _convert_params_units(
            groups[group], expected_units_dim, base_units
        )
    dtype_group_map = {
        "PeriodicTorsionPotential": {
            "container": hoomd.md.dihedral.Periodic,  # Should this be periodic, ask Josh
            "parser": _parse_periodic,
        },
        "OPLSTorsionPotential": {
            "container": hoomd.md.dihedral.OPLS,
            "parser": _parse_opls,
        },
        "RyckaertBellemansTorsionPotential": {
            "container": hoomd.md.dihedral.OPLS,  # RBTorsion will converted to OPLS
            "parser": _parse_rb,
        },
    }
    hoomd_version = hoomd.version.version.split(",")
    if int(hoomd_version[1]) >= 8:
        dtype_group_map["PeriodicTorsionPotential"] = (
            {
                "container": hoomd.md.dihedral.Periodic,
                "parser": _parse_periodic,
            },
        )
    else:
        # Should this be periodic, deprecated starting from 3.8.0
        dtype_group_map["PeriodicTorsionPotential"] = (
            {
                "container": hoomd.md.dihedral.Harmonic,
                "parser": _parse_periodic,
            },
        )

    dihedral_forces = list()
    for group in groups:
        dihedral_forces.append(
            dtype_group_map[group]["parser"](
                container=dtype_group_map[group]["container"](),
                dtypes=groups[group],
            )
        )

    return dihedral_forces


def _parse_improper_forces(top, potential_types, potential_refs, base_units):
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

    def _parse_harmonic(container, itypes):
        for itype in itypes:
            member_types = sort_member_types(itype)
            container.params["-".join(member_types)] = {
                "k": itype.parameters["k"],
                "chi0": itype.parameters["phi_eq"],  # diff nomenclature?
            }
        return container

    unique_dtypes = top.improper_types(
        filter_by=PotentialFilters.UNIQUE_NAME_CLASS
    )
    groups = dict()
    for itype in unique_dtypes:
        group = potential_types[itype]
        if group not in groups:
            groups[group] = [itype]
        else:
            groups[group].append(itype)

    for group in groups:
        expected_units_dim = potential_refs[group][
            "expected_parameters_dimensions"
        ]
        groups[group] = _convert_params_units(
            groups[group], expected_units_dim, base_units
        )
    hoomd_version = hoomd.version.version.split(",")
    if int(hoomd_version[1]) >= 8:
        itype_group_map = {
            "HarmonicImproperPotenial": {
                "container": hoomd.md.dihedral.Periodic,
                "parser": _parse_harmonic,
            },
        }
    else:
        # Should this be periodic, deprecated starting from 3.8.0
        itype_group_map = {
            "HarmonicImproperPotenial": {
                "container": hoomd.md.dihedral.Harmonic,
                "parser": _parse_harmonic,
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


def _validate_base_units(base_units, top):
    """Validate the provided base units, infer units (based on top's positions and masses) if none is provided."""
    ref = {
        "energy": u.dimensions.energy,
        "length": u.dimensions.length,
        "mass": u.dimensions.mass,
    }

    unit_systems = {"MD": MD_UNITS, "AKMA": AKMA_UNITS}
    if isinstance(base_units, str):
        base_units = unit_systems[base_units]
    elif isinstance(base_units, dict):
        for key in base_units:
            if key not in ["energy", "mass", "length"]:
                warnings.warn(
                    "Only base unit will be used during the conversion"
                    "i.e., energy, mass, and length, other units provided"
                    "will not be considered."
                )
            else:
                msg = "{key} is in wrong unit dimension"
                assert base_units[key].dimensions == ref[key], msg

        missing = list()
        for base in ["energy", "mass", "length"]:
            if base not in base_units:
                missing.append(base)
        if missing:
            raise (f"base_units is not fully provided, missing {missing}")
    else:
        base_units = _infer_units(top)
        # base_units = {"length": None, "mass": None, "energy": None}

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


def _convert_params_units(potentials, expected_units_dim, base_units):
    """Convert parameters' units in the potential to that specified in the base_units."""
    converted_potentials = list()
    for potential in potentials:
        converted_params = dict()
        for parameter in potential.parameters:
            unit_dim = expected_units_dim[parameter]
            ind_units = re.sub("[^a-zA-Z]+", " ", unit_dim).split()
            for unit in ind_units:
                unit_dim = unit_dim.replace(unit, str(base_units[unit]))
            converted_params[parameter] = potential.parameters[parameter].to(
                unit_dim
            )
        potential.parameters = converted_params
        converted_potentials.append(potential)
    return converted_potentials


def _generate_pairs_list(top):
    """Return a list of pairs that have non-zero scaling factor."""
    nb_scalings, coulombic_scalings = top.scaling_factors

    pairs12 = list()
    if nb_scalings[0] or coulombic_scalings[0]:
        for bond in top.bonds:
            pairs12.append(sorted(bond.connection_members))

    pairs13 = list()
    if nb_scalings[1] or coulombic_scalings[1]:
        for angle in top.angles:
            pairs13.append(
                sorted(
                    [angle.connection_members[0], angle.connection_members[2]]
                )
            )

    pairs14 = list()
    if nb_scalings[2] or coulombic_scalings[2]:
        for dihedral in top.dihedrals:
            pairs14.append(
                sorted(
                    [
                        dihedral.connection_members[0],
                        dihedral.connection_members[-1],
                    ]
                )
            )
    return pairs12, pairs13, pairs14
