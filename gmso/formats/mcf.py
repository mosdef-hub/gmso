"""Write Cassandra Monte Carlo MCF files."""

import datetime
import warnings

import networkx as nx
import numpy as np
import symengine
import unyt as u

from gmso import __version__
from gmso.core.topology import Topology
from gmso.core.views import PotentialFilters
from gmso.exceptions import GMSOError
from gmso.formats.formats_registry import saves_as
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.utils.compatibility import check_compatibility
from gmso.utils.conversions import (
    convert_opls_to_ryckaert,
    convert_ryckaert_to_fourier,
)

__all__ = ["write_mcf"]

potential_templates = PotentialTemplateLibrary()


@saves_as(".mcf")
def write_mcf(top, filename):
    """Generate a Cassandra MCF from a gmso.core.Topology object.

    The MCF file stores the topology information for a single
    species (i.e., molecule) in the Cassandra Monte Carlo software
    (https://cassandra.nd.edu). The gmso.Topology object provided to
    this function should therefore be for a single molecule with the
    relevant forcefield parameters. One MCF file will be required
    for each unique species in the system.

    Parameters
    ----------
    top : gmso.core.Topology
        Topology object
    filename : str or list
        Path of the output file(s). If a string is provided
        and the Topology object has more than one subtopology,
        an integer suffix will be appended to the string.
        If a list is provided and there is more than one subtopology,
        the names of the output files will be
        each element in the list. The number of element in the list
        should match the
        number of unique subtopologies.

    Notes
    -----
    Atom indexing begins at 1. See
    https://cassandra.nd.edu/index.php/documentation for a complete
    description of the MCF format.
    """
    subtops = []
    for molecule in top.unique_site_labels(name_only=True):
        subtops.append(top.create_subtop("molecule", (molecule, 0)))

    if len(subtops) > 1:
        if len(filename) != len(subtops):
            raise ValueError(
                "write_mcf: Number of filenames must match number of unique species in the Topology object"
            )

    for idx, subtop in enumerate(subtops):
        _check_compatibility(subtop)

        # Identify atoms in rings and Cassandra 'fragments'
        in_ring, frag_list, frag_conn = _id_rings_fragments(subtop)

        if len(subtops) > 1:
            if isinstance(filename, str):
                filename = filename[:-4] + f"_{idx}.mcf"
            elif isinstance(filename, list):
                filename = filename[idx]

        # Now we write the MCF file
        with open(filename, "w") as mcf:
            header = (
                "!***************************************"
                "****************************************\n"
                "!Molecular connectivity file\n"
                "!***************************************"
                "****************************************\n"
                f"!File {filename} written by gmso {__version__} "
                f"at {str(datetime.datetime.now())}\n\n"
            )

            mcf.write(header)
            _write_atom_information(mcf, subtop, in_ring)
            _write_bond_information(mcf, subtop)
            _write_angle_information(mcf, subtop)
            _write_dihedral_information(mcf, subtop)
            _write_improper_information(mcf, subtop)
            _write_fragment_information(mcf, subtop, frag_list, frag_conn)
            _write_intrascaling_information(mcf, subtop)

            # That's all, folks!
            mcf.write("\n\nEND\n")


def _id_rings_fragments(top):
    """Identify the rings and fragments of the molecule.

    Parameters
    ----------
    top : gmso.core.Topology
        Topology object

    Returns
    -------
    in_ring : list
        True for each atom in a ring
    frag_list : list
        Atom ids belonging to each fragment
    frag_conn : list
        Fragment ids of connected fragments
    """
    # Identify atoms in rings
    bond_graph = nx.Graph()
    bond_graph.add_edges_from(
        [
            [
                top.get_index(bond.connection_members[0]),
                top.get_index(bond.connection_members[1]),
            ]
            for bond in top.bonds
        ]
    )
    if len(top.bonds) == 0:
        warnings.warn(
            "No bonds found. Cassandra will interpet this as a rigid species"
        )
        in_ring = [False] * len(top.sites)
        frag_list = []
        frag_conn = []
        return in_ring, frag_list, frag_conn

    # Check if entire molecule is connected. Warn if not.
    if not nx.is_connected(bond_graph):
        raise ValueError(
            "Not all components of the molecule are connected. "
            "MCF files are for a single molecule and thus "
            "everything should be connected through bonds."
        )

    all_rings = nx.cycle_basis(bond_graph)
    in_ring = [False] * bond_graph.number_of_nodes()
    adj_to_ring = [False] * bond_graph.number_of_nodes()
    for ring in all_rings:
        for idx in ring:
            in_ring[idx] = True

    # Identify fragments
    # See Shah and Maginn, JCP, 135, 134121, 2011, doi:10.1063/1.3644939
    frag_list = []
    frag_conn = []

    # First create a neighbor list for each atom
    neigh_dict = {
        i: list(bond_graph.neighbors(i))
        for i in range(bond_graph.number_of_nodes())
    }

    # Handle fused/adjoining rings
    rings_changed = True
    while rings_changed:
        rings_changed = False
        for ring1 in all_rings:
            if rings_changed:
                break
            for ring2 in all_rings:
                if ring1 == ring2:
                    continue
                if len(set(ring1) & set(ring2)) > 0:
                    all_rings.remove(ring1)
                    all_rings.remove(ring2)
                    all_rings.append(list(set(ring1 + ring2)))
                    rings_changed = True
                    break

    # ID fragments which contain a ring
    for ring in all_rings:
        adjacent_atoms = []
        for idx in ring:
            if len(neigh_dict[idx]) > 2:
                adjacent_atoms.append(list(set(neigh_dict[idx]) - set(ring)))
        tmp = filter(None, adjacent_atoms)
        adjacent_atoms = [x for sublist in tmp for x in sublist]
        frag_list.append(ring + adjacent_atoms)
        for idx in adjacent_atoms:
            adj_to_ring[idx] = True
    # Now ID the other fragments
    for idx in neigh_dict:
        if len(neigh_dict[idx]) > 1:
            if in_ring[idx] == True:
                continue
            else:
                frag_list.append([idx] + neigh_dict[idx])
    # Now find connectivity (shared bonds)
    for i in range(len(frag_list)):
        frag1 = frag_list[i]
        for j in range(i + 1, len(frag_list)):
            frag2 = frag_list[j]
            shared_atoms = list(set(frag1) & set(frag2))
            if len(shared_atoms) == 2:
                frag_conn.append([i, j])
            elif len(shared_atoms) > 2:
                warnings.warn(
                    "Fragments share more than two atoms..."
                    "something may be going awry unless there are"
                    "fused rings in your system. See below for details."
                )
                print("Fragment 1 atoms:")
                print(frag1)
                print("Fragment 2 atoms:")
                print(frag2)

    return in_ring, frag_list, frag_conn


def _write_atom_information(mcf, top, in_ring):
    """Write the atoms in the system.

    Parameters
    ----------
    mcf : file object
        The file object of the Cassandra mcf being written
    top : Topology
        Topology object
    in_ring : list
        Boolean for each atom idx True if atom belongs to a ring
    """
    # Based upon Cassandra; updated following 1.2.2 release
    max_element_length = 6
    max_atomtype_length = 20

    sites, names, atypes_list = zip(
        *[(site, site.name, site.atom_type.name) for site in top.sites]
    )

    # Check constraints on atom type length and element name length
    n_unique_names = len(set(names))
    for name in names:
        if len(name) > max_element_length:
            warnings.warn(
                f"Name: {name} will be shortened to {max_element_length}"
                "characters. Please confirm your final MCF."
            )

    # Confirm that shortening names to two characters does not
    # cause two previously unique atom names to become identical.
    names = [name[:max_element_length] for name in names]
    if len(set(names)) < n_unique_names:
        warnings.warn(
            "The number of unique names has been reduced due"
            f"to shortening the name to {max_element_length} characters."
        )

    pfilter = PotentialFilters.UNIQUE_SORTED_NAMES
    n_unique_types = top.atom_types(filter_by=pfilter)

    for type_ in n_unique_types:
        if len(type_.name) > max_atomtype_length:
            warnings.warn(
                f"Type name: {type_.name} will be shortened to "
                f"{max_atomtype_length} characters as "
                f"{type[-max_atomtype_length:]}. Please confirm your final MCF."
            )
    atypes_list = [itype[-max_atomtype_length:] for itype in atypes_list]

    if len(set(atypes_list)) < len(n_unique_types):
        warnings.warn(
            "The number of unique atomtypes has been reduced due to "
            f"shortening the atomtype name to {max_atomtype_length} characters."
        )

    # Check charge neutrality
    net_q = 0.0
    for idx, site in enumerate(range(top.n_sites)):
        net_q += top.sites[idx].charge.in_units(u.elementary_charge).value

    if not np.isclose(net_q, 0.0):
        raise ValueError(
            "Net charge of the system is not zero. "
            "Cassandra MFC requires a neutral system."
        )

    # Detect VDW style
    vdw_styles = set()
    for site in sites:
        vdw_styles.add(_get_vdw_style(site))
    if len(vdw_styles) > 1:
        raise GMSOError(
            "More than one vdw_style detected. "
            "Cassandra only supports MCF files with a single "
            "vdw_style"
        )
    if False in vdw_styles:
        raise GMSOError("Unsupported vdw_style detected.")

    vdw_style = vdw_styles.pop()

    header = (
        "!Atom Format\n"
        "!index type element mass charge vdw_type parameters\n"
        '!vdw_type="LJ", parms=epsilon sigma\n'
        '!vdw_type="Mie", parms=epsilon sigma '
        "repulsion_exponent dispersion_exponent\n"
        "\n# Atom_Info\n"
    )

    mcf.write(header)
    mcf.write("{:d}\n".format(len(top.sites)))
    for idx, site in enumerate(top.sites):
        mcf.write(
            "{:<4d}  "
            "{:<6s}  "
            "{:<2s}  "
            "{:8.4f}  "
            "{:12.8f}  ".format(
                idx + 1,
                atypes_list[idx],
                names[idx],
                site.atom_type.mass.in_units(u.amu).value,
                site.charge.in_units(u.elementary_charge).value,
            )
        )
        if vdw_style == "LJ":
            mcf.write(
                "{:3s}  "
                "{:10.5f}  "
                "{:10.5f}".format(
                    vdw_style,
                    (site.atom_type.parameters["epsilon"] / u.kb)
                    .in_units("K")
                    .value,
                    site.atom_type.parameters["sigma"]
                    .in_units("Angstrom")
                    .value,
                )
            )
        elif vdw_style == "Mie":
            mcf.write(
                "{:3s}  "
                "{:10.5f}  "
                "{:10.5f}  "
                "{:8.3f}  "
                "{:8.3f}".format(
                    vdw_style,
                    (site.atom_type.parameters["epsilon"] / u.kb)
                    .in_units("K")
                    .value,
                    site.atom_type.parameters["sigma"]
                    .in_units("Angstrom")
                    .value,
                    site.atom_type.parameters["n"].value,
                    site.atom_type.parameters["m"].value,
                )
            )
        if in_ring[idx] == True:
            mcf.write("  ring")
        mcf.write("\n")


def _write_bond_information(mcf, top):
    """Write the bonds in the system.

    Parameters
    ----------
    mcf : file object
        The file object of the Cassandra mcf being written
    top : Topology
        Topology object

    """
    mcf.write("\n!Bond Format\n")
    mcf.write(
        "!index i j type parameters\n" + '!type="fixed", parms=bondLength\n'
    )
    mcf.write("\n# Bond_Info\n")
    mcf.write("{:d}\n".format(len(top.bonds)))
    for idx, bond in enumerate(top.bonds):
        mcf.write(
            "{:<4d}  "
            "{:<4d}  "
            "{:<4d}  "
            "{:s}  "
            "{:10.5f}\n".format(
                idx + 1,
                top.get_index(bond.connection_members[0]) + 1,
                top.get_index(bond.connection_members[1]) + 1,
                "fixed",
                bond.connection_type.parameters["r_eq"]
                .in_units(u.Angstrom)
                .value,
            )
        )


def _write_angle_information(mcf, top):
    """Write the angles in the system.

    Parameters
    ----------
    mcf : file object
        The file object of the Cassandra mcf being written
    top : Topology
        Topology object
    """
    header = (
        "\n!Angle Format\n"
        "!index i j k type parameters\n"
        '!type="fixed", parms=equilibrium_angle\n'
        '!type="harmonic", parms=force_constant equilibrium_angle\n'
        "\n# Angle_Info\n"
    )

    mcf.write(header)
    mcf.write("{:d}\n".format(len(top.angles)))
    for idx, angle in enumerate(top.angles):
        mcf.write(
            f"{idx + 1:<4d}  "
            f"{top.get_index(angle.connection_members[0]) + 1:<4d}  "
            f"{top.get_index(angle.connection_members[1]) + 1:<4d}  "
            f"{top.get_index(angle.connection_members[2]) + 1:<4d}  "
        )
        angle_style = _get_angle_style(angle)
        if angle_style == "fixed":
            mcf.write(
                "{:s}  "
                "{:10.5f}\n".format(
                    angle_style,
                    angle.connection_type.parameters["theta_eq"]
                    .in_units(u.degree)
                    .value,
                )
            )
        elif angle_style == "harmonic":
            mcf.write(
                "{:s}  "
                "{:10.5f} "
                "{:10.5f}\n".format(
                    angle_style,
                    (0.5 * angle.connection_type.parameters["k"] / u.kb)
                    .in_units("K/rad**2")
                    .value,  # TODO: k vs. k/2. conversion
                    angle.connection_type.parameters["theta_eq"]
                    .in_units(u.degree)
                    .value,
                )
            )
        else:
            raise GMSOError("Unsupported angle style for Cassandra MCF writer")


def _write_dihedral_information(mcf, top):
    """Write the dihedrals in the system.

    Parameters
    ----------
    mcf : file object
        The file object of the Cassandra mcf being written
    top : Topology
        Topology object
    """
    # Dihedral info
    header = (
        "\n!Dihedral Format\n"
        "!index i j k l type parameters\n"
        '!type="none"\n'
        '!type="CHARMM", parms=a0 a1 delta\n'
        '!type="OPLS", parms=c0 c1 c2 c3\n'
        '!type="harmonic", parms=force_constant equilibrium_dihedral\n'
        "\n# Dihedral_Info\n"
    )

    mcf.write(header)

    mcf.write("{:d}\n".format(len(top.dihedrals)))
    for idx, dihedral in enumerate(top.dihedrals):
        mcf.write(
            "{:<4d}  "
            "{:<4d}  "
            "{:<4d}  "
            "{:<4d}  "
            "{:<4d}  ".format(
                idx + 1,
                top.get_index(dihedral.connection_members[0]) + 1,
                top.get_index(dihedral.connection_members[1]) + 1,
                top.get_index(dihedral.connection_members[2]) + 1,
                top.get_index(dihedral.connection_members[3]) + 1,
            )
        )
        dihedral_style = _get_dihedral_style(dihedral)
        dihedral_style = dihedral_style.upper()
        # If ryckaert, convert to opls
        if dihedral_style == "RYCKAERT":
            dihedral.connection_type = convert_ryckaert_to_fourier(
                dihedral.connection_type
            )
            # The GMSO Fourier style is named OPLS in Cassandra. See
            # https://cassandra-mc.readthedocs.io/en/latest/guides/forcefield.html?highlight=opls#dihedrals

            dihedral_style = "FOURIER"

        if dihedral_style == "OPLS":
            dihedral.connection_type = convert_opls_to_ryckaert(
                dihedral.connection_type
            )
            dihedral.connection_type = convert_ryckaert_to_fourier(
                dihedral.connection_type
            )
            dihedral_style = "FOURIER"

        if dihedral_style == "FOURIER":
            # Cassandra supports the OPLS (GMSO's Fourier) functional form up 4 terms, namely
            # E = a0 + a1 * (1 + cos(phi)) + a2 * (1 - cos(2*phi)) + a3 * (1 + cos(3*phi))
            # The GMSO Fourier potential has terms up to 0.5 * k4 * ( 1 - cos ( 4 * phi))
            # So we need to exclude the last term in the GMSO topology.
            dihedral_style = "OPLS"
            mcf.write(
                "{:s}  "
                "{:10.5f}  "
                "{:10.5f}  "
                "{:10.5f}  "
                "{:10.5f}\n".format(
                    dihedral_style,
                    0.5
                    * dihedral.connection_type.parameters["k0"]
                    .in_units("kJ/mol")
                    .value,
                    0.5
                    * dihedral.connection_type.parameters["k1"]
                    .in_units("kJ/mol")
                    .value,
                    0.5
                    * dihedral.connection_type.parameters["k2"]
                    .in_units("kJ/mol")
                    .value,
                    0.5
                    * dihedral.connection_type.parameters["k3"]
                    .in_units("kJ/mol")
                    .value,
                )
            )
        elif dihedral_style == "CHARMM":
            mcf.write(
                "{:s}  "
                "{:10.5f}  "
                "{:10.5f}  "
                "{:10.5f}\n".format(
                    dihedral_style,
                    dihedral.connection_type.parameters["k"]
                    .in_units("kJ/mol")
                    .value,
                    dihedral.connection_type.parameters["n"],
                    dihedral.connection_type.parameters["phi_eq"]
                    .in_units(u.degrees)
                    .value,
                )
            )
        elif dihedral_style == "HARMONIC":
            mcf.write(
                "{:s}  "
                "{:10.5f}  "
                "{:10.5f}\n".format(
                    dihedral_style.lower(),
                    0.5
                    * dihedral.connection_type.parameters["k"]
                    .in_units("kJ/mol")
                    .value,
                    dihedral.connection_type.parameters["phi_eq"]
                    .in_units(u.degrees)
                    .value,
                )
            )

        else:
            raise GMSOError(
                "Unsupported dihedral style for Cassandra MCF writer"
            )


def _write_improper_information(mcf, top):
    """Write the impropers in the system.

    Parameters
    ----------
    mcf : file object
        The file object of the Cassandra mcf being written
    top : Topology
        Topology object
    """
    header = (
        "\n!Improper Format\n"
        "!index i j k l type parameters\n"
        '!type="harmonic", parms=force_constant equilibrium_improper\n'
        "\n# Improper_Info\n"
    )

    mcf.write(header)
    mcf.write("{:d}\n".format(len(top.impropers)))

    improper_style = "harmonic"
    for i, improper in enumerate(top.impropers):
        mcf.write(
            "{:<4d}  {:<4d}  {:<4d}  {:<4d}  {:<4d}"
            "  {:s}  {:10.5f}  {:10.5f}\n".format(
                i + 1,
                top.get_index(improper.connection_members[0]) + 1,
                top.get_index(improper.connection_members[1]) + 1,
                top.get_index(improper.connection_members[2]) + 1,
                top.get_index(improper.connection_members[3]) + 1,
                improper_style,
                0.5 * improper.connection_type.parameters["k"],
                improper.connection_type.parameters["phi_eq"],
            )
        )


def _write_fragment_information(mcf, top, frag_list, frag_conn):
    """Write the fragments in the molecule.

    Parameters
    ----------
    mcf : file object
        The file object of the Cassandra mcf being written
    top : gmso.core.Topology
        Topology object
    frag_list : list
        Atom ids belonging to each fragment
    frag_conn : list
        Fragment ids of connected fragments

    """
    header = (
        "\n!Fragment Format\n"
        "!index number_of_atoms_in_fragment branch_point other_atoms\n"
        "\n# Fragment_Info\n"
    )

    mcf.write(header)

    # Special cases first
    if len(frag_list) == 0:
        if len(top.sites) == 1:
            mcf.write("1\n")
            mcf.write("1 1 1\n")
        elif len(top.sites) == 2:
            mcf.write("1\n")
            mcf.write("1 2 1 2\n")
        else:
            warnings.warn(
                "More than two atoms present but no fragments identified."
            )
            mcf.write("0\n")
    else:
        mcf.write("{:d}\n".format(len(frag_list)))
        for i, frag in enumerate(frag_list):
            mcf.write("{:d}    {:d}".format(i + 1, len(frag)))
            for idx in frag:
                mcf.write("    {:d}".format(idx + 1))
            mcf.write("\n")

    mcf.write("\n\n# Fragment_Connectivity\n")
    mcf.write("{:d}\n".format(len(frag_conn)))
    for i, conn in enumerate(frag_conn):
        mcf.write(
            "{:d}    {:d}    {:d}\n".format(i + 1, conn[0] + 1, conn[1] + 1)
        )


def _write_intrascaling_information(mcf, top):
    """Write the intramolecular scaling in the molecule.

    Parameters
    ----------
    mcf : file object
        The file object of the Cassandra mcf being written
    lj14 : float
        The 1-4 scaling parameter for LJ interactions
    coul14 : float
        The 1-4 scaling parameter for Coulombic interactions

    """
    nbonded_sf = top.get_lj_scale()
    electstatic_sf = top.get_electrostatics_scale()
    header = (
        "\n!Intra Scaling\n"
        "!vdw_scaling    1-2 1-3 1-4 1-N\n"
        "!charge_scaling 1-2 1-3 1-4 1-N\n"
        "\n# Intra_Scaling\n"
    )

    mcf.write(header)
    mcf.write("{:.4f} {:.4f} {:.4f} 1.0000\n".format(*nbonded_sf))
    mcf.write("{:.4f} {:.4f} {:.4f} 1.0000\n".format(*electstatic_sf))


def _check_compatibility(top):
    """Check Topology object for compatibility with Cassandra MCF format."""
    if not isinstance(top, Topology):
        raise GMSOError("MCF writer requires a Topology object.")
    if not all([site.atom_type for site in top.sites]):
        raise GMSOError(
            "MCF writing not supported without parameterized forcefield."
        )
    accepted_potentials = (
        potential_templates["LennardJonesPotential"],
        potential_templates["MiePotential"],
        potential_templates["FixedBondPotential"],
        potential_templates["HarmonicBondPotential"],
        potential_templates["HarmonicAnglePotential"],
        potential_templates["FixedAnglePotential"],
        potential_templates["PeriodicTorsionPotential"],
        potential_templates["OPLSTorsionPotential"],
        potential_templates["FourierTorsionPotential"],
        potential_templates["RyckaertBellemansTorsionPotential"],
    )
    check_compatibility(top, accepted_potentials)


def _get_vdw_style(site):
    """Return the vdw style."""
    vdw_styles = {
        "LJ": potential_templates["LennardJonesPotential"],
        "Mie": potential_templates["MiePotential"],
    }

    return _get_potential_style(vdw_styles, site.atom_type)


def _get_angle_style(angle):
    """Return the angle style."""
    angle_styles = {
        "harmonic": potential_templates["HarmonicAnglePotential"],
        "fixed": potential_templates["FixedAnglePotential"],
    }

    return _get_potential_style(angle_styles, angle.connection_type)


def _get_dihedral_style(dihedral):
    """Return the dihedral style."""
    dihedral_styles = {
        "charmm": potential_templates["PeriodicTorsionPotential"],
        "harmonic": potential_templates["HarmonicTorsionPotential"],
        "opls": potential_templates["OPLSTorsionPotential"],
        "fourier": potential_templates["FourierTorsionPotential"],
        "ryckaert": potential_templates["RyckaertBellemansTorsionPotential"],
    }

    return _get_potential_style(dihedral_styles, dihedral.connection_type)


def _get_potential_style(styles, potential):
    """Return the potential style."""
    for style, ref in styles.items():
        if ref.independent_variables == potential.independent_variables:
            if symengine.expand(ref.expression - potential.expression) == 0:
                return style
    return False
