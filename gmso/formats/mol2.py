"""Convert to and from a TRIPOS mol2 file."""
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
    The position of atom in the mol2 file is assumed to be in Angstrom.
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
    with open(filename, "r") as f:
        fcontents = f.readlines()
    sections = dict()
    for line in fcontents:
        if "@<TRIPOS>" in line:
            section_key = line.strip("\n")
            sections[section_key] = list()
        else:
            sections[section_key].append(line)

    _parse_site = _parse_lj if site_type == "lj" else _parse_atom
    supported_rti = {
        "@<TRIPOS>ATOM": _parse_site,
        "@<TRIPOS>BOND": _parse_bond,
        "@<TRIPOS>CRYSIN": _parse_box,
        "@<TRIPOS>FF_PBC": _parse_box,
    }
    for section in sections:
        if section not in supported_rti:
            warnings.warn(
                f"The record type indicator {section} is not supported."
            )
        else:
            supported_rti[section](topology, sections[section])

    topology.update_topology()
    # TODO: read in parameters to correct attribute as well. This can be saved in various rti sections.
    return topology


def _parse_lj(top, section):
    """Parse atom of lj style from mol2 file."""
    for line in section:
        cont = line.split()
        position = [float(x) for x in cont[2:5]] * u.Å

        try:
            charge = float(cont[8])
        except IndexError:
            warnings.warn(
                f"No charge was detected for site {cont[1]} with index {cont[0]}"
            )

        atom = Atom(
            name=cont[1],
            position=position.to("nm"),
            charge=charge,
            residue_name=cont[7],
            residue_number=int(cont[6]),
        )
        top.add_site(atom)


def _parse_atom(top, section):
    """Parse atom information from the mol2 file."""
    parse_ele = (
        lambda ele: element_by_symbol(ele)
        if element_by_symbol(ele)
        else element_by_name(ele),
    )
    for line in section:
        cont = line.split()
        position = [float(x) for x in cont[2:5]] * u.Å
        element = parse_ele(cont[5])

        if not element:
            warnings.warn(f"Could not parse the element of {cont[5]}")

        try:
            charge = float(cont[8])
        except IndexError:
            warnings.warn(
                f"No charge was detected for site {cont[1]} with index {cont[0]}"
            )

        atom = Atom(
            name=cont[1],
            position=position.to("nm"),
            element=element,
            charge=charge,
            residue_name=cont[7],
            residue_number=int(cont[6]),
        )
        top.add_site(atom)


def _parse_bond(top, section):
    """Parse bond information from the mol2 file."""
    for line in section:
        cont = line.split()
        bond = Bond(
            connection_members=(
                top.sites[int(cont[1]) - 1],
                top.sites[int(cont[2]) - 1],
            )
        )
        top.add_connection(bond)


def _parse_box(top, section):
    """Parse box information from the mol2 file."""
    if top.box:
        warnings.warn("Topology already has a box")

    for line in section:
        cont = line.split()
        top.box = Box(
            lengths=[float(x) for x in cont[0:3]] * u.Å,
            angles=[float(x) for x in cont[3:6]] * u.degree,
        )
