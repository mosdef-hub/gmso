"""Write a GROMACS topology (.TOP) file."""
import datetime

import unyt as u

from gmso.core.element import element_by_atom_type
from gmso.exceptions import GMSOError
from gmso.formats.formats_registry import saves_as
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.utils.compatibility import check_compatibility


@saves_as(".top")
def write_top(top, filename, top_vars=None):
    """Write a gmso.core.Topology object to a GROMACS topology (.TOP) file."""
    pot_types = _validate_compatibility(top)
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
            "gen-pairs\t"
            "fudgeLJ\t"
            "fudgeQQ\n"
        )
        out_file.write(
            "{0}\t\t\t"
            "{1}\t\t\t"
            "{2}\t\t\t"
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
            "at.num\t\t"
            "mass\t\t"
            "charge\t\t"
            "ptype\t\t"
            "sigma\t\t"
            "epsilon\n"
        )
        for atom_type in top.atom_types:
            out_file.write(
                "{0}\t\t\t"
                "{1}\t\t\t"
                "{2:.5f}\t\t"
                "{3:.5f}\t\t"
                "{4}\t\t\t"
                "{5:.5f}\t\t\t"
                "{6:.5f}\n".format(
                    atom_type.name,
                    _lookup_atomic_number(atom_type),
                    atom_type.mass.in_units(u.amu).value,
                    atom_type.charge.in_units(u.elementary_charge).value,
                    "A",
                    atom_type.parameters["sigma"].in_units(u.nanometer).value,
                    atom_type.parameters["epsilon"]
                    .in_units(u.Unit("kJ/mol"))
                    .value,
                )
            )

        out_file.write("\n[ moleculetype ]\n" "; name\t\tnrexcl\n")

        # TODO: Better parsing of subtops into residues/molecules
        n_unique_subtops = len(set([s.name for s in top.subtops]))
        if n_unique_subtops > 1:
            raise NotImplementedError
        # Treat top without subtops as one residue-like "molecule"
        elif n_unique_subtops == 0:
            out_file.write(
                "{0}\t\t\t"
                "{1}\n\n".format(
                    top.name,
                    top_vars["nrexcl"],  # Typically exclude 3 nearest neighbors
                )
            )
        # TODO: Lookup and join nrexcl from each subtop object
        elif n_unique_subtops == 1:
            out_file.write("{0}\t\t\t" "{1}\n\n".format(top.name, 3))

        out_file.write(
            "[ atoms ]\n"
            "; nr\t\ttype\tresnr\tresidue\t\tatom\tcgnr\tcharge\t\tmass\n"
        )
        for site in top.sites:
            out_file.write(
                "{0}\t\t\t"
                "{1}\t\t"
                "{2}\t\t"
                "{3}\t"
                "{4}\t\t"
                "{5}\t\t"
                "{6:.5f}\t\t"
                "{7:.5f}\n".format(
                    top.get_index(site) + 1,
                    site.atom_type.name,
                    1,  # TODO: subtop idx
                    top.name,  # TODO: subtop.name
                    _lookup_element_symbol(site.atom_type),
                    1,  # TODO: care about charge groups
                    site.charge.in_units(u.elementary_charge).value,
                    site.atom_type.mass.in_units(u.amu).value,
                )
            )

        out_file.write("\n[ bonds ]\n" ";   ai     aj  funct   c0      c1\n")
        for bond in top.bonds:
            out_file.write(
                _write_connection(top, bond, pot_types[bond.connection_type])
            )

        out_file.write(
            "\n[ angles ]\n" ";   ai     aj      ak      funct   c0      c1\n"
        )
        for angle in top.angles:
            out_file.write(
                _write_connection(top, angle, pot_types[angle.connection_type])
            )

        out_file.write(
            "\n[ dihedrals ]\n"
            ";   ai     aj      ak      al  funct   c0      c1      c2\n"
        )
        for dihedral in top.dihedrals:
            out_file.write(
                _write_connection(
                    top, dihedral, pot_types[dihedral.connection_type]
                )
            )

        out_file.write("\n[ system ]\n" "; name\n" "{0}\n\n".format(top.name))

        if len(set([s.name for s in top.subtops])) > 1:
            raise NotImplementedError

        # TODO: Write out atom types for each unique `subtop` in `atoms` section
        # and write out number of molecules in `molecules` section
        # if len(top.subtops) == 0:
        out_file.write(
            "[ molecules ]\n"
            "; molecule\tnmols\n"
            "{0}\t\t{1}".format(top.name, 1)
        )
        # elif len(top.subtops) > 0:
        #    out_file.write(
        #        '[ molecules ]\n'
        #        '; molecule\tnmols\n'
        #        '{0}\t\t{1}'.format(top.subtops[0].name, top.n_subtops)
        #    )


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


def _validate_compatibility(top):
    """Check compatability of topology object with GROMACS TOP format."""
    pot_types = check_compatibility(top, _accepted_potentials())
    return pot_types


def _get_top_vars(top, top_vars):
    """Generate a dictionary of values for the defaults directive."""
    combining_rule_to_gmx = {"lorentz": 2, "geometric": 3}
    default_top_vars = dict()
    default_top_vars["nbfunc"] = 1
    default_top_vars["comb-rule"] = combining_rule_to_gmx[top.combining_rule]
    default_top_vars["gen-pairs"] = "no"
    default_top_vars["fudgeLJ"] = 1
    default_top_vars["fudgeQQ"] = 1
    default_top_vars["nrexcl"] = 3

    if isinstance(top_vars, dict):
        default_top_vars.update(top_vars)

    return default_top_vars


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


def _write_connection(top, connection, potential_name):
    """Worker function to write various connection information."""
    worker_functions = {
        "HarmonicBondPotential": _harmonic_bond_potential_writer,
        "HarmonicAnglePotential": _harmonic_angle_potential_writer,
        "RyckaertBellemansTorsionPotential": _ryckaert_bellemans_torsion_writer,
        "PeriodicTorsionPotential": _periodic_torsion_writer,
    }

    return worker_functions[potential_name](top, connection)


def _harmonic_bond_potential_writer(top, bond):
    """Write harmonic bond information."""
    line = "\t{0}\t{1}\t{2}\t{3:.5f}\t{4:.5f}\n".format(
        top.get_index(bond.connection_members[0]) + 1,
        top.get_index(bond.connection_members[1]) + 1,
        "1",
        bond.connection_type.parameters["r_eq"].in_units(u.nm).value,
        bond.connection_type.parameters["k"]
        .in_units(u.Unit("kJ / (mol*nm**2)"))
        .value,
    )
    return line


def _harmonic_angle_potential_writer(top, angle):
    """Write harmonic angle information."""
    line = "\t{0}\t{1}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\n".format(
        top.get_index(angle.connection_members[0]) + 1,
        top.get_index(angle.connection_members[1]) + 1,
        top.get_index(angle.connection_members[2]) + 1,
        "1",
        angle.connection_type.parameters["theta_eq"].in_units(u.degree).value,
        angle.connection_type.parameters["k"]
        .in_units(u.Unit("kJ/(mol*rad**2)"))
        .value,
    )
    return line


def _ryckaert_bellemans_torsion_writer(top, dihedral):
    """Write Ryckaert-Bellemans Torsion information."""
    line = "\t{0}\t{1}\t{2}\t{3}\t{4}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10:.5f}\n".format(
        top.get_index(dihedral.connection_members[0]) + 1,
        top.get_index(dihedral.connection_members[1]) + 1,
        top.get_index(dihedral.connection_members[2]) + 1,
        top.get_index(dihedral.connection_members[3]) + 1,
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


def _periodic_torsion_writer(top, dihedral):
    """Write periodic torsion information."""
    line = "\t{0}\t{1}\t{2}\t{3}\t{4}\t{5:.5f}\t{6:.5f}\t{7}\n".format(
        top.get_index(dihedral.connection_members[0]) + 1,
        top.get_index(dihedral.connection_members[1]) + 1,
        top.get_index(dihedral.connection_members[2]) + 1,
        top.get_index(dihedral.connection_members[3]) + 1,
        "1",
        dihedral.connection_type.parameters["phi_eq"].in_units(u.degree).value,
        dihedral.connection_type.parameters["k"]
        .in_units(u.Unit("kJ/(mol)"))
        .value,
        dihedral.connection_type.parameters["n"].value,
    )
    return line
