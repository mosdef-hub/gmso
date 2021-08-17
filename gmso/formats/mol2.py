"""Convert to and from a TRIPOS mol2 file."""
# TODO add sources of mol2 files
import os
import warnings
from pathlib import Path

import unyt as u

from gmso import Atom, Bond, Box, Topology
from gmso.core.element import element_by_name, element_by_symbol


def from_mol2(
    filename, site_type="atom"
):  # TODO add flags for information to return
    """Read in a TRIPOS mol2 file format into a gmso topology object.

    Creates a Topology from a mol2 file structure. This will read in the
    topological structure (sites, bonds, and box) information into gmso.
    Note that parameterized information can be found in these objects, but
    will not be converted to the Topology.

    Parameters
    ----------
    filename : string
        path to the file where the mol2 file is stored.
    site_type : string, optional, options:(default:'atom','lj')
        tells the reader to consider the elements saved in the mol2 file, and
        if the type is 'lj', to not try to identify the element of the site,
        instead saving the site name.

    Returns
    -------
    top : gmso.Topology

    Notes
    -----
    It may be common to want to create an mBuild compound from a mol2 file. This is possible
    by installing [mBuild](https://mbuild.mosdef.org/en/stable/index.html)
    and converting using the following python code:

        >>> from gmso.formats.mol2 import from_mol2
        >>> from gmso.external.convert_mbuild import to_mbuild
        >>> top = from_mol2('myfile.mol2')
        >>> mbuild_compound = to_mbuild(top)
    """
    msg = "Provided path to file that does not exist"
    path = Path(filename)
    if not path.exists():
        raise OSError(msg)
    # Initialize topology
    topology = Topology(name=path.name)
    # save the name from the filename
    f = open(path, "r")
    line = f.readline()
    while f:
        # check for header character in line
        if line.startswith("@<TRIPOS>"):
            # if header character in line, send to a function that will direct it properly
            line = parse_record_type_indicator(
                f, line, topology, site_type
            )
        elif line == "":
            break
        else:
            # else, skip to next line
            line = f.readline()
    f.close()
    topology.update_topology()
    # TODO: read in parameters to correct attribute as well. This can be saved in various rti sections.
    return topology


def load_top_sites(f, topology, site_type="atom"):
    """Take a mol2 file section with the heading '<TRIPOS>ATOM' and save to the topology.sites attribute."""
    while True:
        line = f.readline()
        if "@" not in line and not line == "\n":
            line = line.split()
            position = [float(x) for x in line[2:5]] * u.Å
            # TODO: make sure charges are also saved as a unyt value
            # TODO: add validation for element names
            if site_type == "lj":
                element = None
            elif element_by_symbol(line[5]):
                element = element_by_symbol(line[5])
            elif element_by_name(line[5]):
                element = element_by_name(line[5])
            else:
                warnings.warn(
                    "No element detected for site {} with index{}, consider manually adding the element to the topology".format(
                        line[1], len(topology.sites) + 1
                    )
                )
                element = None
            try:
                charge = float(line[8])
            except IndexError:
                warnings.warn(
                    "No charges were detected for site {} with index {}".format(
                        line[1], line[0]
                    )
                )
                charge = None
            atom = Atom(
                name=line[1],
                position=position.to("nm"),
                charge=charge,
                element=element,
            )
            topology.add_site(atom)
        else:
            break
    return line


def load_top_bonds(f, topology, **kwargs):
    """Take a mol2 file section with the heading '@<TRIPOS>BOND' and save to the topology.bonds attribute."""
    while True:
        line = f.readline()
        if "@" not in line and not line == "\n":
            line = line.split()
            bond = Bond(
                connection_members=(
                    topology.sites[int(line[1]) - 1],
                    topology.sites[int(line[2]) - 1],
                )
            )
            topology.add_connection(bond)
        else:
            break
    return line


def load_top_box(f, topology, **kwargs):
    """Take a mol2 file section with the heading '@<TRIPOS>FF_PBC' or '@<TRIPOS>CRYSIN' and save to topology.box."""
    if topology.box:
        warnings.warn(
            "This mol2 file has two boxes to be read in, only reading in one with dimensions {}".format(
                topology.box
            )
        )
        line = f.readline()
        return line
    while True:
        line = f.readline()
        if "@" not in line and not line == "\n":
            line = line.split()
            # TODO: write to box information
            topology.box = Box(
                lengths=[float(x) for x in line[0:3]] * u.Å,
                angles=[float(x) for x in line[3:6]] * u.degree,
            )
        else:
            break
    return line


def parse_record_type_indicator(f, line, topology, site_type):
    """Take a specific record type indicator from a mol2 file format and save to the proper attribute of a gmso topology.

    Supported record type indicators include Atom, Bond, FF_PBC, and CRYSIN.
    """
    supported_rti = {
        "@<TRIPOS>ATOM\n": load_top_sites,
        "@<TRIPOS>BOND\n": load_top_bonds,
        "@<TRIPOS>CRYSIN\n": load_top_box,
        "@<TRIPOS>FF_PBC\n": load_top_box,
    }
    # read in to atom attribute
    try:
        line = supported_rti[line](f, topology, site_type=site_type)
    except KeyError:
        warnings.warn(
            "The record type indicator {} is not supported".format(line)
        )
        line = f.readline()
    return line
