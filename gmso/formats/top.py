"""Write a GROMACS topology (.TOP) file."""
import datetime
from enum import unique

import unyt as u

from gmso.core.dihedral import Dihedral
from gmso.core.element import element_by_atom_type
from gmso.core.improper import Improper
from gmso.core.views import PotentialFilters
from gmso.exceptions import GMSOError
from gmso.formats.formats_registry import saves_as
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.parameterization.molecule_utils import (
    molecule_angles,
    molecule_bonds,
    molecule_dihedrals,
    molecule_impropers,
)
from gmso.utils.compatibility import check_compatibility


@saves_as(".top")
def write_top(top, filename, top_vars=None, simplify_check=False):
    """Write a gmso.core.Topology object to a GROMACS topology (.TOP) file."""
    pot_types = _validate_compatibility(top, simplify_check)
    top_vars = _get_top_vars(top, top_vars)

    with open(filename, "w") as out_file:
        out_file.write(
            "; File {} written by GMSO at {}\n\n".format(
                top.name if top.name is not None else "",
                str(datetime.datetime.now()),
            )
        )
        out_file.write(
            "[ defaults ]\n"
            "; nbfunc\t"
            "comb-rule\t"
            "gen-pairs\t\t"
            "fudgeLJ\t"
            "fudgeQQ\n"
        )
        out_file.write(
            "{0}\t\t\t"
            "{1}\t\t\t"
            "{2}\t\t"
            "{3}\t\t"
            "{4}\n\n".format(
                top_vars["nbfunc"],
                top_vars["comb-rule"],
                top_vars["gen-pairs"],
                top_vars["fudgeLJ"],
                top_vars["fudgeQQ"],
            )
        )

        out_file.write(
            "[ atomtypes ]\n"
            "; name\t\t"
            "at.num\t"
            "mass\t\t"
            "charge\t\t"
            "ptype\t"
            "sigma\t"
            "epsilon\n"
        )

        for atom_type in top.atom_types(PotentialFilters.UNIQUE_NAME_CLASS):
            out_file.write(
                "{0:12s}"
                "{1:4s}"
                "{2:12.5f}"
                "{3:12.5f}\t"
                "{4:4s}"
                "{5:12.5f}"
                "{6:12.5f}\n".format(
                    atom_type.name,
                    str(_lookup_atomic_number(atom_type)),
                    atom_type.mass.in_units(u.amu).value,
                    atom_type.charge.in_units(u.elementary_charge).value,
                    "A",
                    atom_type.parameters["sigma"].in_units(u.nanometer).value,
                    atom_type.parameters["epsilon"]
                    .in_units(u.Unit("kJ/mol"))
                    .value,
                )
            )

        # Define unique molecule by name only
        unique_molecules = _get_unique_molecules(top)

        # Section headers
        headers = {
            "bonds": "\n[ bonds ]\n" ";\tai\t\taj\t\tfunct\t\tc0\t\tc1\n",
            "angles": "\n[ angles ]\n"
            ";\ta\t\taj\t\tak\t\tfunct\t\tc0\t\tc1\n",
            "angle_restraints": (
                "\n[ angle_restraints ]\n"
                ";\tai\t\taj\t\tai\t\tak\t\tfunct\t\ttheta_eq\t\tk\t\tmultiplicity\n"
            ),
            "dihedrals": {
                "RyckaertBellemansTorsionPotential": "\n[ dihedrals ]\n"
                ";\tai\t\taj\t\tak \tal\t\tfunct\t\tc0\t\tc1\t\tc2\t\tc3\t\tc4\n",
                "PeriodicTorsionPotential": "\n[ dihedrals ]\n"
                ";\tai\t\taj\t\tak\t\tal\t\tfunct\t\tphi\t\tk_phi\t\tmulitplicity\n",
            },
            "dihedral_restraints": "\n[ dihedral_restraints ]\n"
            "#ifdef DIHRES\n"
            ";\tai\t\taj\t\tak\t\tal\t\tfunct\t\ttheta_eq\t\tdelta_theta\t\tkd\n",
        }
        for tag in unique_molecules:
            """Write out nrexcl for each unique molecule."""
            out_file.write("\n[ moleculetype ]\n" "; name\tnrexcl\n")

            # TODO: Lookup and join nrexcl from each molecule object
            out_file.write("{0}\t" "{1}\n".format(tag, top_vars["nrexcl"]))

            """Write out atoms for each unique molecule."""
            out_file.write(
                "[ atoms ]\n"
                "; nr\ttype\tresnr\tresidue\t\tatom\tcgnr\tcharge\tmass\n"
            )
            # Each unique molecule need to be reindexed (restarting from 0)
            # The shifted_idx_map is needed to make sure all the atom index used in
            # latter connection sections are acurate
            shifted_idx_map = dict()
            for idx, site in enumerate(unique_molecules[tag]["sites"]):
                shifted_idx_map[top.get_index(site)] = idx
                out_file.write(
                    "\t{0:8s}"
                    "{1:12s}"
                    "{2:8s}"
                    "{3:12s}"
                    "{4:8s}"
                    "{5:4s}"
                    "{6:12.5f}"
                    "{7:12.5f}\n".format(
                        str(idx + 1),
                        site.atom_type.name,
                        str(site.molecule.number + 1 if site.molecule else 1),
                        tag,
                        site.element.symbol if site.element else site.name,
                        "1",  # TODO: care about charge groups
                        site.charge.in_units(u.elementary_charge).value,
                        site.atom_type.mass.in_units(u.amu).value,
                    )
                )

            for conn_group in [
                "bonds",
                "angles",
                "angle_restraints",
                "dihedrals",
                "dihedral_restraints",
                "impropers",
            ]:
                if unique_molecules[tag][conn_group]:
                    if conn_group in ["dihedrals", "impropers"]:
                        proper_groups = {
                            "RyckaertBellemansTorsionPotential": list(),
                            "PeriodicTorsionPotential": list(),
                        }
                        for dihedral in unique_molecules[tag][conn_group]:
                            ptype = pot_types[dihedral.connection_type]
                            proper_groups[ptype].append(dihedral)

                        # Improper use same header as dihedral periodic header
                        if proper_groups["RyckaertBellemansTorsionPotential"]:
                            out_file.write(
                                headers["dihedrals"][
                                    "RyckaertBellemansTorsionPotential"
                                ]
                            )
                            for conn in proper_groups[
                                "RyckaertBellemansTorsionPotential"
                            ]:
                                for line in _write_connection(
                                    top,
                                    conn,
                                    pot_types[conn.connection_type],
                                    shifted_idx_map,
                                ):
                                    out_file.write(line)
                        if proper_groups["PeriodicTorsionPotential"]:
                            out_file.write(
                                headers["dihedrals"]["PeriodicTorsionPotential"]
                            )
                            for conn in proper_groups[
                                "PeriodicTorsionPotential"
                            ]:
                                for line in _write_connection(
                                    top,
                                    conn,
                                    pot_types[conn.connection_type],
                                    shifted_idx_map,
                                ):
                                    out_file.write(line)
                    else:
                        out_file.write(headers[conn_group])
                        for conn in unique_molecules[tag][conn_group]:
                            out_file.write(
                                _write_connection(
                                    top,
                                    conn,
                                    pot_types[conn.connection_type],
                                    shifted_idx_map,
                                )
                            )
                    if conn_group == "dihedral_restraints":
                        out_file.write("#endif DIHRES\n")

        out_file.write("\n[ system ]\n" "; name\n" "{0}\n\n".format(top.name))

        out_file.write("[ molecules ]\n" "; molecule\tnmols\n")
        for tag in unique_molecules:
            out_file.write(
                "{0}\t{1}\n".format(tag, len(unique_molecules[tag]["subtags"]))
            )


def _accepted_potentials():
    """List of accepted potentials that GROMACS can support."""
    templates = PotentialTemplateLibrary()
    lennard_jones_potential = templates["LennardJonesPotential"]
    harmonic_bond_potential = templates["HarmonicBondPotential"]
    harmonic_angle_potential = templates["HarmonicAnglePotential"]
    periodic_torsion_potential = templates["PeriodicTorsionPotential"]
    rb_torsion_potential = templates["RyckaertBellemansTorsionPotential"]
    accepted_potentials = [
        lennard_jones_potential,
        harmonic_bond_potential,
        harmonic_angle_potential,
        periodic_torsion_potential,
        rb_torsion_potential,
    ]
    return accepted_potentials


def _validate_compatibility(top, simplify_check):
    """Check compatability of topology object with GROMACS TOP format."""
    pot_types = check_compatibility(top, _accepted_potentials(), simplify_check)
    return pot_types


def _get_top_vars(top, top_vars):
    """Generate a dictionary of values for the defaults directive."""
    combining_rule_to_gmx = {"lorentz": 2, "geometric": 3}
    default_top_vars = dict()
    default_top_vars["nbfunc"] = 1  # modify this to check for lj or buckingham
    default_top_vars["comb-rule"] = combining_rule_to_gmx[top.combining_rule]
    default_top_vars["gen-pairs"] = "no"
    default_top_vars["fudgeLJ"] = top.scaling_factors[0][2]
    default_top_vars["fudgeQQ"] = top.scaling_factors[1][2]
    default_top_vars["nrexcl"] = 3

    if isinstance(top_vars, dict):
        default_top_vars.update(top_vars)

    return default_top_vars


def _get_unique_molecules(top):
    unique_molecules = {
        tag: {
            "subtags": list(),
        }
        for tag in top.unique_site_labels("molecule", name_only=True)
    }

    for molecule in top.unique_site_labels("molecule", name_only=False):
        unique_molecules[molecule.name]["subtags"].append(molecule)

    if len(unique_molecules) == 0:
        unique_molecules[top.name] = dict()
        unique_molecules[top.name]["subtags"] = [top.name]
        unique_molecules[top.name]["sites"] = list(top.sites)
        unique_molecules[top.name]["bonds"] = list(top.bonds)
        unique_molecules[top.name]["angles"] = list(top.angles)
        unique_molecules[top.name]["angle_restraints"] = list(
            angle for angle in top.angles if angle.restraint
        )
        unique_molecules[top.name]["dihedrals"] = list(top.angles)
        unique_molecules[top.name]["dihedral_restraints"] = list(
            dihedral for dihedral in top.dihedrals if dihedral.restraint
        )
        unique_molecules[molecule.name]["impropers"] = list(top.impropers)

    else:
        for tag in unique_molecules:
            molecule = unique_molecules[tag]["subtags"][0]
            unique_molecules[tag]["sites"] = list(
                top.iter_sites(key="molecule", value=molecule)
            )
            unique_molecules[tag]["bonds"] = list(molecule_bonds(top, molecule))
            unique_molecules[tag]["angles"] = list(
                molecule_angles(top, molecule)
            )
            unique_molecules[tag]["angle_restraints"] = list(
                angle
                for angle in molecule_angles(top, molecule)
                if angle.restraint
            )
            unique_molecules[tag]["dihedrals"] = list(
                molecule_dihedrals(top, molecule)
            )
            unique_molecules[tag]["dihedral_restraints"] = list(
                dihedral
                for dihedral in molecule_dihedrals(top, molecule)
                if dihedral.restraint
            )
            unique_molecules[tag]["impropers"] = list(
                molecule_impropers(top, molecule)
            )
    return unique_molecules


def _lookup_atomic_number(atom_type):
    """Look up an atomic_number based on atom type information, 0 if non-element type."""
    try:
        element = element_by_atom_type(atom_type)
        return element.atomic_number
    except GMSOError:
        return 0


def _lookup_element_symbol(atom_type):
    """Look up an atomic_number based on atom type information, 0 if non-element type."""
    try:
        element = element_by_atom_type(atom_type)
        return element.symbol
    except GMSOError:
        return "X"


def _write_connection(top, connection, potential_name, shifted_idx_map):
    """Worker function to write various connection information."""
    worker_functions = {
        "HarmonicBondPotential": _harmonic_bond_potential_writer,
        "HarmonicAnglePotential": _harmonic_angle_potential_writer,
        "RyckaertBellemansTorsionPotential": _ryckaert_bellemans_torsion_writer,
        "PeriodicTorsionPotential": _periodic_torsion_writer,
    }

    return worker_functions[potential_name](top, connection, shifted_idx_map)


def _harmonic_bond_potential_writer(top, bond, shifted_idx_map):
    """Write harmonic bond information."""
    line = "\t{0:8s}{1:8s}{2:4s}{3:15.5f}{4:15.5f}\n".format(
        str(shifted_idx_map[top.get_index(bond.connection_members[0])] + 1),
        str(shifted_idx_map[top.get_index(bond.connection_members[1])] + 1),
        "1",
        bond.connection_type.parameters["r_eq"].in_units(u.nm).value,
        bond.connection_type.parameters["k"]
        .in_units(u.Unit("kJ / (mol*nm**2)"))
        .value,
    )
    return line


def _harmonic_angle_potential_writer(top, angle, shifted_idx_map):
    """Write harmonic angle information."""
    line = "\t{0:8s}{1:8s}{2:8s}{3:4s}{4:15.5f}{5:15.5f}\n".format(
        str(shifted_idx_map[top.get_index(angle.connection_members[0])] + 1),
        str(shifted_idx_map[top.get_index(angle.connection_members[1])] + 1),
        str(shifted_idx_map[top.get_index(angle.connection_members[2])] + 1),
        "1",
        angle.connection_type.parameters["theta_eq"].in_units(u.degree).value,
        angle.connection_type.parameters["k"]
        .in_units(u.Unit("kJ/(mol*rad**2)"))
        .value,
    )
    return line


def _ryckaert_bellemans_torsion_writer(top, dihedral, shifted_idx_map):
    """Write Ryckaert-Bellemans Torsion information."""
    line = "\t{0:8s}{1:8s}{2:8s}{3:8s}{4:4s}{5:15.5f}{6:15.5f}{7:15.5f}{8:15.5f}{9:15.5f}{10:15.5f}\n".format(
        str(shifted_idx_map[top.get_index(dihedral.connection_members[0])] + 1),
        str(shifted_idx_map[top.get_index(dihedral.connection_members[1])] + 1),
        str(shifted_idx_map[top.get_index(dihedral.connection_members[2])] + 1),
        str(shifted_idx_map[top.get_index(dihedral.connection_members[3])] + 1),
        "3",
        dihedral.connection_type.parameters["c0"]
        .in_units(u.Unit("kJ/mol"))
        .value,
        dihedral.connection_type.parameters["c1"]
        .in_units(u.Unit("kJ/mol"))
        .value,
        dihedral.connection_type.parameters["c2"]
        .in_units(u.Unit("kJ/mol"))
        .value,
        dihedral.connection_type.parameters["c3"]
        .in_units(u.Unit("kJ/mol"))
        .value,
        dihedral.connection_type.parameters["c4"]
        .in_units(u.Unit("kJ/mol"))
        .value,
        dihedral.connection_type.parameters["c5"]
        .in_units(u.Unit("kJ/mol"))
        .value,
    )
    return line


def _periodic_torsion_writer(top, dihedral, shifted_idx_map):
    """Write periodic torsion information."""
    if isinstance(dihedral, Dihedral):
        if len(dihedral.connection_type.parameters["phi_eq"]) == 1:
            # Normal dihedral
            layers, funct = 1, "1"
        else:
            # Layered/Multiple dihedral
            layers, funct = (
                len(dihedral.connection_type.parameters["phi_eq"]),
                "9",
            )
    elif isinstance(dihedral, Improper):
        layers, funct = 1, "4"
    else:
        raise TypeError(f"Type {type(dihedral)} not supported.")

    lines = list()
    for i in range(layers):
        line = "\t{0:8s}{1:8s}{2:8s}{3:8s}{4:4s}{5:15.5f}{6:15.5f}{7:4.5f}\n".format(
            str(
                shifted_idx_map[top.get_index(dihedral.connection_members[0])]
                + 1
            ),
            str(
                shifted_idx_map[top.get_index(dihedral.connection_members[1])]
                + 1
            ),
            str(
                shifted_idx_map[top.get_index(dihedral.connection_members[2])]
                + 1
            ),
            str(
                shifted_idx_map[top.get_index(dihedral.connection_members[3])]
                + 1
            ),
            funct,
            dihedral.connection_type.parameters["phi_eq"][i]
            .in_units(u.degree)
            .value,
            dihedral.connection_type.parameters["k"][i]
            .in_units(u.Unit("kJ/(mol)"))
            .value,
            dihedral.connection_type.parameters["n"][i].value,
        )
        lines.append(line)
    return lines


def _write_restraint(top, connection, type, shifted_idx_map):
    """Worker function to write various connection restraint information."""
    worker_functions = {
        "angle": _angle_restraint_writer,
        "dihedral": _dihedral_restraint_writer,
    }

    return worker_functions[type](top, connection, shifted_idx_map)


def _angle_restraint_writer(top, angle, shifted_idx_map):
    """Write angle restraint information."""
    line = (
        "\t{0:8s}{1:8s}{2:8s}{3:8s}{4:4s}{5:15.5f}{6:15.5f}{7:4.5f}\n".format(
            str(
                shifted_idx_map[top.get_index(angle.connection_members[1])] + 1
            ),
            str(
                shifted_idx_map[top.get_index(angle.connection_members[0])] + 1
            ),
            str(
                shifted_idx_map[top.get_index(angle.connection_members[1])] + 1
            ),
            str(
                shifted_idx_map[top.get_index(angle.connection_members[2])] + 1
            ),
            "1",
            angle.restraint["theta_eq"].in_units(u.degree).value,
            angle.restraint["k"].in_units(u.Unit("kJ/(mol)")).value,
            angle.restraint["n"],
        )
    )
    return line


def _dihedral_restraint_writer(top, dihedral, shifted_idx_map):
    """Write dihedral restraint information."""
    line = (
        "\t{0:8s}{1:8s}{2:8s}{3:8s}{4:4s}{5:15.5f}{6:15.5f}{7:4.5f}\n".format(
            str(
                shifted_idx_map[top.get_index(dihedral.connection_members[0])]
                + 1
            ),
            str(
                shifted_idx_map[top.get_index(dihedral.connection_members[1])]
                + 1
            ),
            str(
                shifted_idx_map[top.get_index(dihedral.connection_members[2])]
                + 1
            ),
            str(
                shifted_idx_map[top.get_index(dihedral.connection_members[3])]
                + 1
            ),
            "1",
            dihedral.restraint["phi_eq"].in_units(u.degree).value,
            dihedral.restraint["delta_phi"].in_units(u.degree).value,
            dihedral.restraint["k"].in_units(u.Unit("kJ/(mol * rad**2)")).value,
        )
    )
    return line
