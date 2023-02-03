"""Read and write LAMMPS data files."""
from __future__ import division

import datetime
import warnings

import numpy as np
import unyt as u
from sympy import simplify, sympify
from unyt.array import allclose_units
from pathlib import Path

from gmso.core.angle import Angle
from gmso.core.angle_type import AngleType
from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.bond import Bond
from gmso.core.bond_type import BondType
from gmso.core.views import PotentialFilters as pfilters
from gmso.core.box import Box
from gmso.core.element import element_by_mass
from gmso.core.topology import Topology
from gmso.formats.formats_registry import loads_as, saves_as
from gmso.lib.potential_templates import PotentialTemplateLibrary
from gmso.utils.conversions import (
    convert_opls_to_ryckaert,
    convert_ryckaert_to_opls,
)
from gmso.utils.compatibility import check_compatibility
from gmso.utils.decorators import mark_WIP


@saves_as(".lammps", ".lammpsdata", ".data")
@mark_WIP("Testing in progress")
def write_lammpsdata(top, filename, atom_style="full", unit_style="real", strict_potentials=False, strict_units=False):
    """Output a LAMMPS data file.

    Outputs a LAMMPS data file in the 'full' atom style format.
    Assumes use of 'real' units.
    See http://lammps.sandia.gov/doc/atom_style.html for more information on atom styles.

    Parameters
    ----------
    Topology : `Topology`
        A Topology Object
    filename : str
        Path of the output file
    atom_style : str, optional, default='full'
        Defines the style of atoms to be saved in a LAMMPS data file.
        The following atom styles are currently supported: 'full', 'atomic', 'charge', 'molecular'
        see http://lammps.sandia.gov/doc/atom_style.html for more information on atom styles.
    strict : bool, optional, default False
        Tells the writer how to treat conversions. If strict=False, then check for conversions
        of unit styles in #TODO

    Notes
    -----
    See http://lammps.sandia.gov/doc/2001/data_format.html for a full description of the LAMMPS data format.
    This is a work in progress, as only atoms, masses, and atom_type information can be written out.

    Some of this function has been adopted from `mdtraj`'s support of the LAMMPSTRJ trajectory format.
    See https://github.com/mdtraj/mdtraj/blob/master/mdtraj/formats/lammpstrj.py for details.

    """
    # TODO: Support atomstyles ["atomic", "charge", "molecular", "full"]
    if atom_style not in ["full"]:
        raise ValueError(
            'Atom style "{}" is invalid or is not currently supported'.format(
                atom_style
            )
        )

    # TODO: Support various unit styles ["metal", "si", "cgs", "electron", "micro", "nano"]
    if unit_style not in ["real"]:
        raise ValueError(
            'Unit style "{}" is invalid or is not currently supported'.format(
                unit_style
            )
        )
    # Use gmso unit packages to get into correct lammps formats
    default_unit_maps = {"real": "TODO"}
    default_parameter_maps = { # Add more as needed
        "dihedrals":"OPLSTorsionPotential",
        "angles":"HarmonicAnglePotential",
        "bonds":"HarmonicBondPotential",
        #"atoms":"LennardJonesPotential",
        #"electrostatics":"CoulombicPotential"
    }

    # TODO: Use strict_x to validate depth of topology checking
    if strict_units:
        _validate_unit_compatibility(top, default_unit_maps[unit_style])
    else:
        top = _try_default_unit_conversions(top, default_unit_maps[unit_style])


    if strict_potentials:
        print("I'm strict about potential forms")
        _validate_potential_compatibility(top)
    else:
        top = _try_default_potential_conversions(top, default_parameter_maps)

    # TODO: improve handling of various filenames
    path = Path(filename)
    if not path.parent.exists():
         msg = "Provided path to file that does not exist"
         raise FileNotFoundError(msg)

    with open(path, "w") as out_file:
        _write_header(out_file, top, atom_style)
        _write_box(out_file, top)
        if top.is_typed(): #TODO: should this be is_fully_typed?
            _write_atomtypes(out_file, top)
            _write_pairtypes(out_file, top)
            if top.bonds: _write_bondtypes(out_file, top)
            if top.angles: _write_angletypes(out_file, top)
            if top.dihedrals: _write_dihedraltypes(out_file, top)
            if top.impropers: _write_impropertypes(out_file, top)

        _write_site_data(out_file, top, atom_style)
        for conn in ["bonds", "angles", "dihedrals", "impropers"]:
            connIter = getattr(top, conn)
            if connIter: _write_conn_data(out_file, top, connIter, conn)


@loads_as(".lammps", ".lammpsdata", ".data")
@mark_WIP("Testing in progress")
def read_lammpsdata(
    filename, atom_style="full", unit_style="real", potential="lj"
):
    """Read in a lammps data file as a GMSO topology.

    Parameters
    ----------
    filename : str
        LAMMPS data file
    atom_style : str, optional, default='full'
        Inferred atom style defined by LAMMPS
    potential: str, optional, default='lj'
        Potential type defined in data file

    Returns
    -------
    top : GMSO Topology
        A GMSO Topology object

    Notes
    -----
    See http://lammps.sandia.gov/doc/2001/data_format.html for a full description of the LAMMPS data format.

    This is a work in progress, as only several atom styles, potential styles, and unit styes are currently supported.

    Currently supporting the following atom styles: 'full'

    Currently supporting the following unit styles: 'real'

    Currently supporting the following potential styles: 'lj'

    Proper dihedrals can be read in but is currently not tested.

    Currently not supporting improper dihedrals.

    """
    # TODO: This whole function probably needs to be revamped
    # TODO: Add argument to ask if user wants to infer bond type
    top = Topology()

    # Validate 'atom_style'
    if atom_style not in ["full"]:
        raise ValueError(
            'Atom Style "{}" is invalid or is not currently supported'.format(
                atom_style
            )
        )

    # Validate 'unit_style'
    if unit_style not in ["real"]:
        raise ValueError(
            'Unit Style "{}" is invalid or is not currently supported'.format(
                unit_style
            )
        )

    # Parse box information
    _get_box_coordinates(filename, unit_style, top)
    # Parse atom type information
    top, type_list = _get_ff_information(filename, unit_style, top)
    # Parse atom information
    _get_atoms(filename, top, unit_style, type_list)
    # Parse connection (bonds, angles, dihedrals) information
    # TODO: Add more atom styles
    if atom_style in ["full"]:
        _get_connection(filename, top, unit_style, connection_type="bond")
        _get_connection(filename, top, unit_style, connection_type="angle")

    top.update_topology()

    return top


def get_units(unit_style):
    """Get units for specific LAMMPS unit style."""
    # Need separate angle units for harmonic force constant and angle
    unit_style_dict = {
        "real": {
            "mass": u.g / u.mol,
            "distance": u.angstrom,
            "energy": u.kcal / u.mol,
            "angle_k": u.radian,
            "angle": u.degree,
            "charge": u.elementary_charge,
        }
    }

    return unit_style_dict[unit_style]


def _get_connection(filename, topology, unit_style, connection_type):
    """Parse connection types."""
    with open(filename, "r") as lammps_file:
        types = False
        for i, line in enumerate(lammps_file):
            if connection_type in line.split():
                n_connection_types = int(line.split()[0])
                types = True
            if connection_type.capitalize() in line.split():
                break
    if types == False:
        return topology
    connection_type_lines = open(filename, "r").readlines()[
        i + 2 : i + n_connection_types + 2
    ]
    connection_type_list = list()
    for line in connection_type_lines:
        if connection_type == "bond":
            c_type = BondType(name=line.split()[0])
            # Multiply 'k' by 2 since LAMMPS includes 1/2 in the term
            c_type.parameters["k"] = (
                float(line.split()[1])
                * u.Unit(
                    get_units(unit_style)["energy"]
                    / get_units(unit_style)["distance"] ** 2
                )
                * 2
            )
            c_type.parameters["r_eq"] = float(line.split()[2]) * (
                get_units(unit_style)["distance"]
            )
        elif connection_type == "angle":
            c_type = AngleType(name=line.split()[0])
            # Multiply 'k' by 2 since LAMMPS includes 1/2 in the term
            c_type.parameters["k"] = (
                float(line.split()[1])
                * u.Unit(
                    get_units(unit_style)["energy"]
                    / get_units(unit_style)["angle_k"] ** 2
                )
                * 2
            )
            c_type.parameters["theta_eq"] = float(line.split()[2]) * u.Unit(
                get_units(unit_style)["angle"]
            )

        connection_type_list.append(c_type)

    with open(filename, "r") as lammps_file:
        for i, line in enumerate(lammps_file):
            if connection_type + "s" in line.split():
                n_connections = int(line.split()[0])
            if connection_type.capitalize() + "s" in line.split():
                break
    connection_lines = open(filename, "r").readlines()[
        i + 2 : i + n_connections + 2
    ]
    # Determine number of sites to generate
    if connection_type == "bond":
        n_sites = 2
    elif connection_type == "angle":
        n_sites = 3
    else:
        n_sites = 4
    for i, line in enumerate(connection_lines):
        site_list = list()
        for j in range(n_sites):
            site = topology.sites[int(line.split()[j + 2]) - 1]
            site_list.append(site)
        if connection_type == "bond":
            connection = Bond(
                connection_members=site_list,
                bond_type=connection_type_list[int(line.split()[1]) - 1],
            )
        elif connection_type == "angle":
            connection = Angle(
                connection_members=site_list,
                angle_type=connection_type_list[int(line.split()[1]) - 1],
            )
        topology.add_connection(connection)

    return topology


def _get_atoms(filename, topology, unit_style, type_list):
    """Parse the atom information in the LAMMPS data file."""
    with open(filename, "r") as lammps_file:
        for i, line in enumerate(lammps_file):
            if "atoms" in line.split():
                n_atoms = int(line.split()[0])
            if "Atoms" in line.split():
                break
    atom_lines = open(filename, "r").readlines()[i + 2 : i + n_atoms + 2]
    for line in atom_lines:
        atom_line = line.split()
        atom_type = atom_line[2]
        charge = u.unyt_quantity(
            float(atom_line[3]), get_units(unit_style)["charge"]
        )
        coord = u.angstrom * u.unyt_array(
            [float(atom_line[4]), float(atom_line[5]), float(atom_line[6])]
        )
        site = Atom(
            charge=charge,
            position=coord,
            atom_type=type_list[int(atom_type) - 1],
        )
        element = element_by_mass(site.atom_type.mass.value)
        site.name = element.name
        site.element = element
        topology.add_site(site)

    return topology


def _get_box_coordinates(filename, unit_style, topology):
    """Parse box information."""
    with open(filename, "r") as lammps_file:
        for line in lammps_file:
            if "xlo" in line.split():
                break
        x_line = line.split()
        y_line = lammps_file.readline().split()
        z_line = lammps_file.readline().split()

        x = float(x_line[1]) - float(x_line[0])
        y = float(y_line[1]) - float(y_line[0])
        z = float(z_line[1]) - float(z_line[0])

        # Check if box is triclinic
        tilts = lammps_file.readline().split()
        if "xy" in tilts:
            xy = float(tilts[0])
            xz = float(tilts[1])
            yz = float(tilts[2])

            xhi = float(x_line[1]) - np.max([0.0, xy, xz, xy + xz])
            xlo = float(x_line[0]) - np.min([0.0, xy, xz, xy + xz])
            yhi = float(y_line[1]) - np.max([0.0, yz])
            ylo = float(y_line[0]) - np.min([0.0, yz])
            zhi = float(z_line[1])
            zlo = float(z_line[0])

            lx = xhi - xlo
            ly = yhi - ylo
            lz = zhi - zlo

            c = np.sqrt(lz**2 + xz**2 + yz**2)
            b = np.sqrt(ly**2 + xy**2)
            a = lx

            alpha = np.arccos((yz * ly + xy * xz) / (b * c))
            beta = np.arccos(xz / c)
            gamma = np.arccos(xy / b)

            # Box Information
            lengths = u.unyt_array([a, b, c], get_units(unit_style)["distance"])
            angles = u.unyt_array([alpha, beta, gamma], u.radian)
            angles.to(get_units(unit_style)["angle"])
            topology.box = Box(lengths, angles)
        else:
            # Box Information
            lengths = u.unyt_array([x, y, z], get_units(unit_style)["distance"])
            topology.box = Box(lengths)

        return topology


def _get_ff_information(filename, unit_style, topology):
    """Parse atom-type information."""
    with open(filename, "r") as lammps_file:
        types = False
        for i, line in enumerate(lammps_file):
            if "atom" in line:
                n_atomtypes = int(line.split()[0])
                types = True
            elif "Masses" in line:
                break
    if types == False:
        return topology
    mass_lines = open(filename, "r").readlines()[i + 2 : i + n_atomtypes + 2]
    type_list = list()
    for line in mass_lines:
        atom_type = AtomType(
            name=line.split()[0],
            mass=float(line.split()[1]) * get_units(unit_style)["mass"],
        )
        type_list.append(atom_type)

    with open(filename, "r") as lammps_file:
        for i, line in enumerate(lammps_file):
            if "Pair" in line:
                break
    # Need to figure out if we're going have mixing rules printed out
    # Currently only reading in LJ params
    pair_lines = open(filename, "r").readlines()[i + 2 : i + n_atomtypes + 2]
    for i, pair in enumerate(pair_lines):
        if len(pair.split()) == 3:
            type_list[i].parameters["sigma"] = (
                float(pair.split()[2]) * get_units(unit_style)["distance"]
            )
            type_list[i].parameters["epsilon"] = (
                float(pair.split()[1]) * get_units(unit_style)["energy"]
            )
        elif len(pair.split()) == 4:
            warnings.warn("Currently not reading in mixing rules")

    return topology, type_list

def _accepted_potentials():
    """List of accepted potentials that LAMMPS can support."""
    templates = PotentialTemplateLibrary()
    lennard_jones_potential = templates["LennardJonesPotential"]
    harmonic_bond_potential = templates["HarmonicBondPotential"]
    harmonic_angle_potential = templates["HarmonicAnglePotential"]
    periodic_torsion_potential = templates["PeriodicTorsionPotential"]
    fourier_torsion_potential = templates["FourierTorsionPotential"]
    accepted_potentialsList = [
        lennard_jones_potential,
        harmonic_bond_potential,
        harmonic_angle_potential,
        periodic_torsion_potential,
        fourier_torsion_potential,
    ]
    return accepted_potentialsList


def _validate_potential_compatibility(top):
    """Check compatability of topology object potentials with LAMMPSDATA format."""
    pot_types = check_compatibility(top, _accepted_potentials())
    return pot_types

def _validate_unit_compatibility(top, unitSet):
    """Check compatability of topology object units with LAMMPSDATA format."""
    # TODO: Check to make sure all units are in the correct format
    return True

# All writer worker function belows
def _write_header(out_file, top, atom_style):
    """Write Lammps file header"""
    out_file.write(
        "{} written by topology at {} using the GMSO LAMMPS Writer\n\n".format(
            top.name if top.name is not None else "",
            str(datetime.datetime.now()),
        )
    )
    out_file.write("{:d} atoms\n".format(top.n_sites))
    if atom_style in ["full", "molecular"]:
        out_file.write("{:d} bonds\n".format(top.n_bonds))
        out_file.write("{:d} angles\n".format(top.n_angles))
        out_file.write("{:d} dihedrals\n".format(top.n_dihedrals))
        out_file.write("{:d} impropers\n".format(top.n_impropers))

    # TODO: allow users to specify filter_by syntax
    out_file.write("\n{:d} atom types\n".format(len(top.atom_types(filter_by=pfilters.UNIQUE_NAME_CLASS))))
    if top.n_bonds > 0:
        out_file.write("{:d} bond types\n".format(len(top.bond_types(filter_by=pfilters.UNIQUE_NAME_CLASS))))
    if top.n_angles > 0:
        out_file.write("{:d} angle types\n".format(len(top.angle_types(filter_by=pfilters.UNIQUE_NAME_CLASS))))
    if top.n_dihedrals > 0:
        out_file.write("{:d} dihedral types\n".format(len(top.dihedral_types(filter_by=pfilters.UNIQUE_NAME_CLASS))))
    if top.n_impropers > 0:
        out_file.write("{:d} improper types\n".format(len(top.improper_types(filter_by=pfilters.UNIQUE_NAME_CLASS))))

    out_file.write("\n")

def _write_box(out_file, top):
    """Write GMSO Topology box to LAMMPS file."""
    # TODO: unit conversions
    if allclose_units(
        top.box.angles,
        u.unyt_array([90, 90, 90], "degree"),
        rtol=1e-5,
        atol=1e-8,
    ):
        top.box.lengths.convert_to_units(u.angstrom)
        for i, dim in enumerate(["x", "y", "z"]):
            out_file.write(
                "{0:.6f} {1:.6f} {2}lo {2}hi\n".format(
                    0, top.box.lengths.value[i], dim
                )
            )
        out_file.write("0.000000 0.000000 0.000000 xy xz yz\n")
    else:
        top.box.lengths.convert_to_units(u.angstrom)
        top.box.angles.convert_to_units(u.radian)
        vectors = top.box.get_vectors()
        a, b, c = top.box.lengths
        alpha, beta, gamma = top.box.angles

        lx = a
        xy = b * np.cos(gamma)
        xz = c * np.cos(beta)
        ly = np.sqrt(b**2 - xy**2)
        yz = (b * c * np.cos(alpha) - xy * xz) / ly
        lz = np.sqrt(c**2 - xz**2 - yz**2)

        xhi = vectors[0][0]
        yhi = vectors[1][1]
        zhi = vectors[2][2]
        xy = vectors[1][0]
        xz = vectors[2][0]
        yz = vectors[2][1]
        xlo = u.unyt_array(0, xy.units)
        ylo = u.unyt_array(0, xy.units)
        zlo = u.unyt_array(0, xy.units)

        xlo_bound = xlo + u.unyt_array(
            np.min([0.0, xy, xz, xy + xz]), xy.units
        )
        xhi_bound = xhi + u.unyt_array(
            np.max([0.0, xy, xz, xy + xz]), xy.units
        )
        ylo_bound = ylo + u.unyt_array(np.min([0.0, yz]), xy.units)
        yhi_bound = yhi + u.unyt_array(np.max([0.0, yz]), xy.units)
        zlo_bound = zlo
        zhi_bound = zhi

        out_file.write(
            "{0:.6f} {1:.6f} xlo xhi\n".format(
                xlo_bound.value, xhi_bound.value
            )
        )
        out_file.write(
            "{0:.6f} {1:.6f} ylo yhi\n".format(
                ylo_bound.value, yhi_bound.value
            )
        )
        out_file.write(
            "{0:.6f} {1:.6f} zlo zhi\n".format(
                zlo_bound.value, zhi_bound.value
            )
        )
        out_file.write(
            "{0:.6f} {1:.6f} {2:.6f} xy xz yz\n".format(
                xy.value, xz.value, yz.value
            )
        )

def _write_atomtypes(out_file, top):
    """Write out atomtypes in GMSO topology to LAMMPS file."""
    # TODO: Get a dictionary of indices and atom types
    # TODO: Allow for unit conversions for the unit styles
    out_file.write("\nMasses\n")
    out_file.write(f"#\tmass ({top.sites[0].mass.units})\n")
    atypesView = top.atom_types(filter_by=pfilters.UNIQUE_NAME_CLASS)
    for atom_type in top.atom_types(filter_by=pfilters.UNIQUE_NAME_CLASS):
        out_file.write(
            "{:d}\t{:.6f}\t# {}\n".format(
                atypesView.index(atom_type) + 1,
                atom_type.mass.in_units(u.g / u.mol).value,
                atom_type.name,
            )
        )

def _write_pairtypes(out_file, top):
    """Write out pair interaction to LAMMPS file."""
    # TODO: Modified cross-interactions
    # TODO: Utilize unit styles and nonbonded equations properly
    # Pair coefficients
    test_atmtype = top.sites[0].atom_type
    out_file.write(f"\nPair Coeffs # lj\n") # TODO: This should be pulled from the test_atmtype
    # TODO: use unit style specified for writer
    param_labels = map(lambda x: f"{x} ({test_atmtype.parameters[x].units})", test_atmtype.parameters)
    out_file.write("#\t" + "\t".join(param_labels) + "\n")
    for idx, param in enumerate(top.atom_types(filter_by=pfilters.UNIQUE_NAME_CLASS)):
        # TODO: grab expression from top
        out_file.write(
            "{}\t{:7.5f}\t\t{:7.5f}\t\t#{}\n".format(
                idx + 1,
                param.parameters["epsilon"]
                .in_units(u.Unit("kcal/mol"))
                .value,
                param.parameters["sigma"]
                .in_units(u.angstrom)
                .value,
                param.name,
            )
        )

def _write_bondtypes(out_file, top):
    """Write out bonds to LAMMPS file."""
    # TODO: Make sure to perform unit conversions
    # TODO: Use any accepted lammps parameters
    test_bontype = top.bonds[0].bond_type
    out_file.write(f"\nBond Coeffs #{test_bontype.name}\n")
    param_labels = map(lambda x: f"{x} ({test_bontype.parameters[x].units})", test_bontype.parameters)
    out_file.write("#\t" + "\t".join(param_labels) + "\n")
    for idx, bond_type in enumerate(top.bond_types(filter_by=pfilters.UNIQUE_NAME_CLASS)):
        out_file.write(
            "{}\t{:7.5f}\t{:7.5f}\t\t# {}\t{}\n".format(
                idx + 1,
                bond_type.parameters["k"]
                .in_units(u.Unit("kcal/mol/angstrom**2"))
                .value,
                bond_type.parameters["r_eq"]
                .in_units(u.Unit("angstrom"))
                .value,
                bond_type.member_types[0],
                bond_type.member_types[1]
            )
        )

def _write_angletypes(out_file, top):
    """Write out angles to LAMMPS file."""
    # TODO: Make sure to perform unit conversions
    # TODO: Use any accepted lammps parameters
    test_angtype = top.angles[0].angle_type
    out_file.write(f"\nAngle Coeffs #{test_angtype.name}\n")
    param_labels = map(lambda x: f"{x} ({test_angtype.parameters[x].units})", test_angtype.parameters)
    out_file.write("#\t" + "\t".join(param_labels) + "\n")
    for idx, angle_type in enumerate(top.angle_types(filter_by=pfilters.UNIQUE_NAME_CLASS)):
        out_file.write(
            "{}\t{:7.5f}\t{:7.5f}\n".format(
                idx + 1,
                angle_type.parameters["k"]
                .in_units(u.Unit("kcal/mol/radian**2"))
                .value,
                angle_type.parameters["theta_eq"]
                .in_units(u.Unit("degree"))
                .value,
            )
        )

def _write_dihedraltypes(out_file, top):
    """Write out dihedrals to LAMMPS file."""
    test_dihtype = top.dihedrals[0].dihedral_type
    print(test_dihtype.parameters)
    out_file.write(f"\nDihedral Coeffs #{test_dihtype.name}\n")
    param_labels = map(lambda x: f"{x} ({test_dihtype.parameters[x].units})", test_dihtype.parameters)
    out_file.write("#\t" + "\t".join(param_labels) + "\n")
    #out_file.write("#\tf1(kcal/mol)\tf2(kcal/mol)\tf3(kcal/mol)\tf4(kcal/mol)\n")
    #out_file.write(f"#\tk ({test_dihtype.parameters[0].units})\t\tthetaeq ({test_dihtype.parameters[1].units}})\n") #check for unit styles
    for idx, dihedral_type in enumerate(top.dihedral_types(filter_by=pfilters.UNIQUE_NAME_CLASS)):
        out_file.write(
            "{}\t{:8.5f}\t{:8.5f}\t{:8.5f}\t{:8.5f}\n".format(
                idx + 1,
                dihedral_type.parameters["k1"]
                .in_units(u.Unit("kcal/mol"))
                .value,
                dihedral_type.parameters["k2"]
                .in_units(u.Unit("kcal/mol"))
                .value,
                dihedral_type.parameters["k3"]
                .in_units(u.Unit("kcal/mol"))
                .value,
                dihedral_type.parameters["k4"]
                .in_units(u.Unit("kcal/mol"))
                .value,
            )
        )

def _write_impropertypes(out_file, top):
    """Write out impropers to LAMMPS file."""
    # TODO: Make sure to perform unit conversions
    # TODO: Use any accepted lammps parameters
    test_imptype = top.impropers[0].improper_type
    out_file.write(f"\nImproper Coeffs #{test_imptype.name}\n")
    param_labels = map(lambda x: f"{x} ({test_imptype.parameters[x].units})", test_imptype.parameters)
    out_file.write("#\t" + "\t".join(param_labels) + "\n")
    for idx, improper_type in enumerate(top.improper_types(filter_by=pfilters.UNIQUE_NAME_CLASS)):
        out_file.write(
            "{}\t{:7.5f}\t{:7.5f}\n".format(
                idx + 1,
                improper_type.parameters["k"]
                .in_units(u.Unit("kcal/mol"))
                .value,
                improper_type.parameters["chieq"]
                .in_units(u.Unit("kcal/mol"))
                .value,
            )
        )

def _write_site_data(out_file, top, atom_style):
    """Write atomic positions and charges to LAMMPS file.."""
    # TODO: Allow for unit system to be passed through
    out_file.write("\nAtoms\n\n")
    if atom_style == "atomic":
        atom_line = "{index:d}\t{type_index:d}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n"
    elif atom_style == "charge":
        atom_line = "{index:d}\t{type_index:d}\t{charge:.6f}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n"
    elif atom_style == "molecular":
        atom_line = "{index:d}\t{zero:d}\t{type_index:d}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n"
    elif atom_style == "full":
        atom_line = "{index:d}\t{zero:d}\t{type_index:d}\t{charge:.6f}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n"

    # TODO: test for speedups in various looping methods
    for i, site in enumerate(top.sites):
        out_file.write(
            atom_line.format(
                index=top.sites.index(site) + 1,
                type_index=top.atom_types(filter_by=pfilters.UNIQUE_NAME_CLASS).equality_index(site.atom_type) + 1,
                zero=0, #What is this zero?
                charge=site.charge.to(u.elementary_charge).value,
                x=site.position[0].in_units(u.angstrom).value,
                y=site.position[1].in_units(u.angstrom).value,
                z=site.position[2].in_units(u.angstrom).value,
            )
        )

def _write_conn_data(out_file, top, connIter, connStr):
    """Write all connections to LAMMPS datafile"""
    # TODO: Test for speedups in various looping methods
    # TODO: Allow for unit system passing
    # TODO: Validate that all connections are written in the correct order
    out_file.write(f"\n{connStr.capitalize()}\n\n")
    indexList = list(map(id, getattr(top, connStr[:-1] + '_types')(filter_by=pfilters.UNIQUE_NAME_CLASS)))
    print(f"Indexed list for {connStr} is {indexList}")
    for i, conn in enumerate(getattr(top, connStr)):
        print(f"{connStr}: id:{id(conn.connection_type)} of form {conn.connection_type}")
        typeStr = f"{i+1:d}\t{getattr(top, connStr[:-1] + '_types')(filter_by=pfilters.UNIQUE_NAME_CLASS).equality_index(conn.connection_type) + 1:1}\t"
        indexStr = "\t".join(map(lambda x: str(top.sites.index(x)+1), conn.connection_members))
        out_file.write(typeStr + indexStr + "\n")

def _try_default_potential_conversions(top, potentialsDict):
    # TODO: Docstrings
    return top.convert_potential_styles(potentialsDict)

def _try_default_unit_conversions(top, unitSet):
    # TODO: Docstrings
    try:
        return top # TODO: Remote this once implemented
        top = top.convert_unit_styles(unitSet)
    except:
        raise ValueError(
             'Unit style "{}" cannot be converted from units used in potential expressions. Check the forcefield for consisten units'.format(
                 unit_style
             )
        )
    return top
