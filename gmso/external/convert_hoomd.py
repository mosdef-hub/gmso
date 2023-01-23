"""Convert GMSO Topology to GSD snapshot."""
from __future__ import division

import itertools
import json
import re
import statistics
import warnings

import hoomd
import numpy as np
import unyt as u
from unyt.array import allclose_units

from gmso.core.bond import Bond
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
from gmso.utils.sorting import natural_sort

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

    if top.n_bonds > 0:
        _parse_bond_information(gsd_snapshot, top)
    if top.n_angles > 0:
        _parse_angle_information(gsd_snapshot, top)
    if top.n_dihedrals > 0:
        _parse_dihedral_information(gsd_snapshot, top)
    if top.n_impropers > 0:
        _parse_improper_information(gsd_snapshot, top)

    return gsd_snapshot


def _parse_particle_information(
    gsd_snapshot,
    top,
    base_units,
    rigid_bodies,
    shift_coords,
    box_lengths,
):
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
    masses = u.unyt_array(masses).to(base_units["mass"])

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
        4.0 * np.pi * e0 * 1 * base_units["length"] * 1 * base_units["energy"]
    ) ** 2

    gsd_snapshot.particles.N = top.n_sites
    gsd_snapshot.particles.position = xyz.in_units(base_units["length"]).value
    gsd_snapshot.particles.types = sorted((set(types)))
    gsd_snapshot.particles.typeid = typeids
    gsd_snapshot.particles.mass = masses.in_units(base_units["mass"]).value
    gsd_snapshot.particles.charge = charges / charge_factor

    if rigid_bodies:
        warnings.warn(
            "Rigid bodies detected, but not yet implemented for GSD",
            NotYetImplementedWarning,
        )


def _parse_bond_information(gsd_snapshot, top):
    """Write the bonds in the system.

    Parameters
    ----------
    gsd_snapshot :
        The file object of the GSD file being written
    top : gmso.Topology
        Topology object holding system information

    """
    gsd_snapshot.bonds.N = top.n_bonds
    warnings.warn(f"{top.n_bonds} bonds detected")
    bond_groups = []
    bond_typeids = []
    bond_types = []

    for bond in top.bonds:
        t1, t2 = list(bond.connection_members)
        if all([t1.atom_type, t2.atom_type]):
            _t1 = t1.atom_type.name
            _t2 = t2.atom_type.name
        else:
            _t1 = t1.name
            _t2 = t2.nam
        _t1, _t2 = sorted([_t1, _t2], key=lambda x: x)
        bond_type = "-".join((_t1, _t2))
        bond_types.append(bond_type)
        bond_groups.append(sorted([top.sites.index(t1), top.sites.index(t2)]))

    unique_bond_types = list(set(bond_types))
    bond_typeids = [unique_bond_types.index(i) for i in bond_types]
    gsd_snapshot.bonds.types = unique_bond_types
    gsd_snapshot.bonds.typeid = bond_typeids
    gsd_snapshot.bonds.group = bond_groups
    warnings.warn(f"{len(unique_bond_types)} unique bond types detected")


def _parse_angle_information(gsd_snapshot, top):
    """Write the angles in the system.

    Parameters
    ----------
    gsd_snapshot :
        The file object of the GSD file being written
    top : gmso.Topology
        Topology object holding system information

    """
    gsd_snapshot.angles.N = top.n_angles
    unique_angle_types = set()
    angle_typeids = []
    angle_groups = []
    angle_types = []

    for angle in top.angles:
        t1, t2, t3 = list(angle.connection_members)
        if all([t1.atom_type, t2.atom_type, t3.atom_type]):
            _t1, _t3 = sorted(
                [t1.atom_type.name, t3.atom_type.name], key=natural_sort
            )
            _t2 = t2.atom_type.name
        else:
            _t1, _t3 = sorted([t1.name, t3.name], key=natural_sort)
            _t2 = t2.name

        angle_type = "-".join((_t1, _t2, _t3))
        angle_types.append(angle_type)
        angle_groups.append(
            (top.sites.index(t1), top.sites.index(t2), top.sites.index(t3))
        )

    unique_angle_types = list(set(angle_types))
    angle_typeids = [unique_angle_types.index(i) for i in angle_types]
    gsd_snapshot.angles.types = unique_angle_types
    gsd_snapshot.angles.typeid = angle_typeids
    gsd_snapshot.angles.group = angle_groups

    warnings.warn(f"{top.n_angles} angles detected")
    warnings.warn(f"{len(unique_angle_types)} unique angle types detected")


def _parse_dihedral_information(gsd_snapshot, top):
    """Write the dihedrals in the system.

    Parameters
    ----------
    gsd_snapshot :
        The file object of the GSD file being written
    top : gmso.Topology
        Topology object holding system information

    """
    gsd_snapshot.dihedrals.N = top.n_dihedrals
    dihedral_groups = []
    dihedral_types = []

    for dihedral in top.dihedrals:
        t1, t2, t3, t4 = list(dihedral.connection_members)
        if all([t.atom_type for t in [t1, t2, t3, t4]]):
            _t1, _t4 = sorted(
                [t1.atom_type.name, t4.atom_type.name], key=natural_sort
            )
            _t3 = t3.atom_type.name
            _t2 = t2.atom_type.name
        else:
            _t1, _t4 = sorted([t1.name, t4.name], key=natural_sort)
            _t2 = t2.name
            _t3 = t3.name

        if [_t2, _t3] == sorted([_t2, _t3], key=natural_sort):
            dihedral_type = "-".join((_t1, _t2, _t3, _t4))
        else:
            dihedral_type = "-".join((_t4, _t3, _t2, _t1))

        dihedral_types.append(dihedral_type)
        dihedral_groups.append(
            (
                top.sites.index(t1),
                top.sites.index(t2),
                top.sites.index(t3),
                top.sites.index(t4),
            )
        )

    unique_dihedral_types = list(set(dihedral_types))
    dihedral_typeids = [unique_dihedral_types.index(i) for i in dihedral_types]
    gsd_snapshot.dihedrals.types = unique_dihedral_types
    gsd_snapshot.dihedrals.typeid = dihedral_typeids
    gsd_snapshot.dihedrals.group = dihedral_groups

    warnings.warn(f"{top.n_dihedrals} dihedrals detected")
    warnings.warn(
        f"{len(unique_dihedral_types)} unique dihedral types detected"
    )


def _parse_improper_information(gsd_snapshot, top):
    """Write the dihedrals in the system.

    Parameters
    ----------
    gsd_snapshot :
        The file object of the GSD file being written
    top : gmso.Topology
        Topology object holding system information

    """
    gsd_snapshot.impropers.N = top.n_impropers
    improper_groups = []
    improper_types = []

    for improper in top.impropers:
        t1, t2, t3, t4 = list(improper.connection_members)
        if all([t.atom_type for t in [t1, t2, t3, t4]]):
            _t1, _t4 = sorted(
                [t1.atom_type.name, t4.atom_type.name], key=natural_sort
            )
            _t3 = t3.atom_type.name
            _t2 = t2.atom_type.name
        else:
            _t1, _t4 = sorted([t1.name, t4.name], key=natural_sort)
            _t2 = t2.name
            _t3 = t3.name

        if [_t2, _t3] == sorted([_t2, _t3], key=natural_sort):
            improper_type = "-".join((_t1, _t2, _t3, _t4))
        else:
            improper_type = "-".join((_t4, _t3, _t2, _t1))

        improper_types.append(improper_type)
        improper_groups.append(
            (
                top.sites.index(t1),
                top.sites.index(t2),
                top.sites.index(t3),
                top.sites.index(t4),
            )
        )

    unique_improper_types = list(set(improper_types))
    improper_typeids = [unique_improper_types.index(i) for i in improper_types]
    gsd_snapshot.impropers.types = unique_improper_types
    gsd_snapshot.impropers.typeid = improper_typeids
    gsd_snapshot.impropers.group = improper_groups

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


def to_hoomd_forcefield(top, nlist_buffer=0.4, base_units=None):
    """Convert the potential portion of a typed GMSO to hoomd forces.

    Parameters
    ----------
    top : gmso.Topology
        The typed topology to be converted
    nlist_buffer : float, optional, default=0.4
        Neighborlist buffer for simulation cell. Its unit is the same as that
        used to defined GMSO Topology Box.
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
        "pairs": _parse_nonbonded_forces(
            top, nlist_buffer, potential_types, potential_refs, base_units
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
    from gmso.utils.compatibility import check_compatibility

    templates = PotentialTemplateLibrary()
    lennard_jones_potential = templates["LennardJonesPotential"]
    harmonic_bond_potential = templates["HarmonicBondPotential"]
    harmonic_angle_potential = templates["HarmonicAnglePotential"]
    periodic_torsion_potential = templates["PeriodicTorsionPotential"]
    opls_torsion_potential = templates["OPLSTorsionPotential"]
    rb_torsion_potential = templates["RyckaertBellemansTorsionPotential"]
    accepted_potentials = [
        lennard_jones_potential,
        harmonic_bond_potential,
        harmonic_angle_potential,
        periodic_torsion_potential,
        opls_torsion_potential,
        rb_torsion_potential,
    ]
    potential_types = check_compatibility(top, accepted_potentials)
    return potential_types


def _parse_nonbonded_forces(
    top, nlist_buffer, potential_types, potential_refs, base_units
):
    """Parse nonbonded forces."""
    # Set up helper methods to parse different nonbonded forces.
    def _parse_lj(container, atypes, combining_rule):
        """Parse LJ forces."""
        for atype1, atype2 in itertools.combinations_with_replacement(
            atypes, 2
        ):
            comb_epsilon = statistics.geometric_mean(
                [atype1.parameters["epsilon"], atype2.parameters["epsilon"]]
            )
            if top.combining_rule == "lorentz":
                comb_sigma = np.mean(
                    [atype1.parameters["sigma"], atype2.parameters["sigma"]]
                )
            elif top.combining_rule == "geometric":
                comb_sigma = statistics.geometric_mean(
                    [atype1.parameters["sigma"], atype2.parameters["sigma"]]
                )
            else:
                raise ValueError(
                    f"Invalid combining rule provided ({combining_rule})"
                )

            container.params[(atype1.name, atype2.name)] = {
                "sigma": comb_sigma,
                "epsilon": comb_epsilon,
            }

        return container

    def _parse_buckingham(container, atypes):
        return None

    def _parse_lj0804(container, atypes):
        return None

    def _parse_lj1208(container, atypes):
        return None

    def _parse_mie(container, atypes):
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
            groups[group], expected_units_dim, base_units
        )

    atype_group_map = {
        "LennardJonesPotential": {
            "container": hoomd.md.pair.LJ,
            "parser": _parse_lj,
        },
        "BuckinghamPotential": {
            "container": hoomd.md.pair.Buckingham,
            "parser": _parse_buckingham,
        },
        "MiePotential": {"container": hoomd.md.pair.Mie, "parser": _parse_mie},
    }

    nlist = hoomd.md.nlist.Cell(exclusions=["bond", "1-3"], buffer=nlist_buffer)

    nbonded_forces = list()
    for group in groups:
        nbonded_forces.append(
            atype_group_map[group]["parser"](
                container=atype_group_map[group]["container"](nlist=nlist),
                atypes=groups[group],
                combining_rule=top.combining_rule,
            )
        )

    return nbonded_forces


def _parse_bond_forces(top, potential_types, potential_refs, base_units):
    """Parse bond forces."""

    def _parse_harmonic(container, btypes):
        for btype in btypes:
            # TODO: Unit conversion
            container.params["-".join(btype.member_types)] = {
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
    """Parse angle forces."""

    def _parse_harmonic(container, agtypes):
        for agtype in agtypes:
            # TODO: Unit conversion
            container.params["-".join(agtype.member_types)] = {
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

    base_units["angle"] = u.degree
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
    """Parse dihedral forces."""

    def _parse_periodic(container, dtypes):
        for dtype in dtypes:
            # TODO: Unit conversion
            container.params["-".join(dtype.member_types)] = {
                "k": dtype.parameters["k"],
                "d": 1,
                "n": dtype.parameters["n"],
                "phi0": dtype.parameters["phi_eq"],
            }
        return container

    def _parse_opls(container, dtypes):
        for dtype in dtypes:
            # TODO: Unit conversion
            # TODO: The range of ks is mismatched (GMSO go from k0 to k5)
            # May need to do a check that k0 == k5 == 0 or raise a warning
            container.params[dtype.member_types] = {
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
            # TODO: Unit conversion
            # TODO: The range of ks is mismatched (GMSO go from k0 to k5)
            # May need to do a check that k0 == k5 == 0 or raise a warning
            container.params["-".join(dtype.member_types)] = {
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
            "container": hoomd.md.dihedral.Harmonic,  # Should this be periodic, ask Josh
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
    """Parse improper forces."""

    def _parse_harmonic(container, itypes):
        for itype in itypes:
            # TODO: Unit conversion
            container.params["-".join(itype.member_types)] = {
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

    itype_group_map = {
        "HarmonicImproperPotenial": {
            "container": hoomd.md.dihedral.Harmonic,  # Should this be periodic, ask Josh
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
    """Try to infer unit from top."""
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
