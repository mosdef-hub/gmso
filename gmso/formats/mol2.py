"""Convert to and from a TRIPOS mol2 file."""
# TODO add sources of mol2 files
import warnings
from pathlib import Path

import unyt as u

from gmso import Atom, Bond, Box, Topology
from gmso.core.element import element_by_name, element_by_symbol
from gmso.formats.formats_registry import loads_as


@loads_as(".mol2")
def from_mol2(filename, site_type="atom"):
    """Read in a TRIPOS mol2 file format into a gmso topology object.

    Creates a Topology from a mol2 file structure. This will read in the
    topological structure (sites, bonds, and box) information into gmso.
    Note that parameterized information can be found in these objects, but
    will not be converted to the Topology.

    Parameters
    ----------
    filename : string
        path to the file where the mol2 file is stored.
    site_type : string ('atom' or 'lj'), default='atom'
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

        >>> from gmso import Topology
        >>> from gmso.external.convert_mbuild import to_mbuild
        >>> top = Topology.load('myfile.mol2')
        >>> mbuild_compound = to_mbuild(top)
    """
    mol2path = Path(filename)
    if not mol2path.exists():
        msg = "Provided path to file that does not exist"
        raise FileNotFoundError(msg)
    # Initialize topology
    topology = Topology(name=mol2path.stem)
    # save the name from the filename
    with open(mol2path, "r") as f:
        line = f.readline()
        while f:
            # check for header character in line
            if line.strip().startswith("@<TRIPOS>"):
                # if header character in line, send to a function that will direct it properly
                line = _parse_record_type_indicator(
                    f, line, topology, site_type
                )
            elif line == "":
                # check for the end of file
                break
            else:
                # else, skip to next line
                line = f.readline()
    topology.update_topology()
    # TODO: read in parameters to correct attribute as well. This can be saved in various rti sections.
    return topology


def _load_top_sites(f, topology, site_type="atom"):
    """Take a mol2 file section with the heading '<TRIPOS>ATOM' and save to the topology.sites attribute.

    Parameters
    ----------
    f : file pointer
        pointer file where the mol2 file is stored. The pointer must be at the head of the rti for that
        `@<TRIPOS>ATOM` section.
    topology : gmso.Topology
        topology to save the site information to.
    site_type : string ('atom' or 'lj'), default='atom'
        tells the reader to consider the elements saved in the mol2 file, and
        if the type is 'lj', to not try to identify the element of the site,
        instead saving the site name.

    Returns
    -------
    line : string
         returns the last line of the `@<TRIPOS>ATOM` section, and this is where the file pointer (`f`)
         will now point to.

    Notes
    -----
    Will modify the topology in place with the relevant site information. Indices will be appended to any
    current site information.

    """
    while True:
        line = f.readline()
        if _is_end_of_rti(line):
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
                residue_name=line[7],
                residue_number=int(line[6]),
            )
            topology.add_site(atom)
        else:
            break
    return line


def _load_top_bonds(f, topology, **kwargs):
    """Take a mol2 file section with the heading '@<TRIPOS>BOND' and save to the topology.bonds attribute."""
    while True:
        line = f.readline()
        if _is_end_of_rti(line):
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


def _load_top_box(f, topology, **kwargs):
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
        if _is_end_of_rti(line):
            line = line.split()
            # TODO: write to box information
            topology.box = Box(
                lengths=[float(x) for x in line[0:3]] * u.Å,
                angles=[float(x) for x in line[3:6]] * u.degree,
            )
        else:
            break
    return line


def _parse_record_type_indicator(f, line, topology, site_type):
    """Take a specific record type indicator (RTI) from a mol2 file format and save to the proper attribute of a gmso topology.

    Supported record type indicators include Atom, Bond, FF_PBC, and CRYSIN.
    """
    supported_rti = {
        "@<TRIPOS>ATOM": _load_top_sites,
        "@<TRIPOS>BOND": _load_top_bonds,
        "@<TRIPOS>CRYSIN": _load_top_box,
        "@<TRIPOS>FF_PBC": _load_top_box,
    }
    # read in to atom attribute
    try:
        line = supported_rti[line.strip()](f, topology, site_type=site_type)
    except KeyError:
        warnings.warn(
            "The record type indicator {} is not supported. Skipping current section and moving to the next RTI header.".format(
                line
            )
        )
        line = f.readline()
    return line


def _is_end_of_rti(line):
    """Check if line in an rti is at the end of the section."""
    return (
        "@" not in line
        and not line == "\n"
        and line
        and not line.strip().startswith("#")
    )
