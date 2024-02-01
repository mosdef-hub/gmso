"""Convert to and from a TRIPOS mol2 file."""

import itertools
import warnings
from pathlib import Path

import unyt as u

from gmso import Atom, Bond, Box, Topology
from gmso.abc.abstract_site import MoleculeType, ResidueType
from gmso.core.element import element_by_name, element_by_symbol
from gmso.formats.formats_registry import loads_as


@loads_as(".mol2")
def read_mol2(filename, site_type="atom", verbose=False):
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
    verbose : bool, optional, default=False
        If True, raise warnings for any assumptions made during the parsing.

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

    sections = {"Meta": list()}
    section_key = "Meta"  # Used to parse the meta info at top of the file
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
        "@<TRIPOS>MOLECULE": _parse_molecule,
    }
    for section in sections:
        if section not in supported_rti:
            warnings.warn(
                f"The record type indicator {section} is not supported. "
                "Skipping current section and moving to the next RTI header."
            )
        else:
            supported_rti[section](topology, sections[section], verbose)

    # TODO: read in parameters to correct attribute as well. This can be saved in various rti sections.
    return topology


def _parse_lj(top, section, verbose):
    """Parse atom of lj style from mol2 file."""
    for line in section:
        if line.strip():
            content = line.split()
            position = [float(x) for x in content[2:5]] * u.Å

            try:
                charge = float(content[8])
            except IndexError:
                if verbose:
                    warnings.warn(
                        f"No charge was detected for site {content[1]} with index {content[0]}"
                    )
                charge = None

            atom = Atom(
                name=content[1],
                position=position.to("nm"),
                charge=charge,
                residue=(content[7], int(content[6])),
            )
            top.add_site(atom)


def _parse_atom(top, section, verbose):
    """Parse atom information from the mol2 file."""

    def parse_ele(*symbols):
        methods = [element_by_name, element_by_symbol]
        elem = None
        for symbol, method in itertools.product(symbols, methods):
            elem = method(symbol)
            if elem:
                break
        return elem

    for line in section:
        if line.strip():
            content = line.split()
            position = [float(x) for x in content[2:5]] * u.Å
            element = parse_ele(content[5], content[1])

            if not element and verbose:
                warnings.warn(
                    f"No element detected for site {content[1]} with index {content[0]}, "
                    "consider manually adding the element to the topology"
                )

            try:
                charge = float(content[8])
            except IndexError:
                if verbose:
                    warnings.warn(
                        f"No charge was detected for site {content[1]} with index {content[0]}"
                    )
                charge = None
            molecule = top.label if top.__dict__.get("label") else top.name
            atom = Atom(
                name=content[1],
                position=position.to("nm"),
                element=element,
                charge=charge,
                residue=ResidueType(content[7], int(content[6])),
                molecule=MoleculeType(molecule, 0),
            )
            top.add_site(atom)


def _parse_bond(top, section, verbose):
    """Parse bond information from the mol2 file."""
    for line in section:
        if line.strip():
            content = line.split()
            bond = Bond(
                connection_members=tuple(
                    sorted(
                        [
                            top.sites[int(content[1]) - 1],
                            top.sites[int(content[2]) - 1],
                        ],
                        key=lambda x: top.get_index(x),
                    )
                )
            )
            top.add_connection(bond)


def _parse_box(top, section, verbose):
    """Parse box information from the mol2 file."""
    if top.box:
        warnings.warn(
            f"This mol2 file has two boxes to be read in, only reading in one with dimensions {top.box}"
        )

    for line in section:
        if line.strip():
            content = line.split()
            top.box = Box(
                lengths=[float(x) for x in content[0:3]] * u.Å,
                angles=[float(x) for x in content[3:6]] * u.degree,
            )


def _parse_molecule(top, section, verbose):
    """Parse molecule information from the mol2 file."""
    top.label = str(section[0].strip())
