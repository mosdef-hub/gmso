"""Convert to and from a TRIPOS mol2 file."""

import datetime
import itertools
import warnings
from pathlib import Path

import unyt as u

import gmso
from gmso import Atom, Bond, Box, Topology
from gmso.abc.abstract_site import Molecule, Residue
from gmso.core.element import element_by_name, element_by_symbol
from gmso.formats.formats_registry import loads_as, saves_as


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


@saves_as(".mol2")
def write_mol2(top, filename, n_decimals=3):
    with open(filename, "w") as out_file:
        out_file.write(
            "{} written by GMSO {} at {}\n".format(
                top.name if top.name is not None else "",
                gmso.__version__,
                str(datetime.datetime.now()),
            )
        )
        _write_molecule_info(top, out_file)
        _write_sites_info(top, out_file)
        _write_bonds_info(top, out_file)
        _write_box_info(top, out_file)


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
                residue=Residue(name=content[7], number=int(content[6])),
                molecule=Molecule(name=molecule, number=0),
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


def _write_sites_info(top, out_file):
    """Write site information to ATOM section."""
    # TODO: Create rules to make sure nothing is too long, so that it cuts off.
    # ATOM top.sites
    # @<TRIPOS>ATOM
    # 1 C          -0.7600    1.1691   -0.0005 C.ar    1  BENZENE       0.000
    out_file.write("@<TRIPOS>ATOM\n")
    for index, site in enumerate(top.sites):
        lineList = [
            str(index + 1),
            site.element.symbol,
            *map(str, site.position.value.round(3)),
            site.atom_type.name if site.atom_type else site.name,
            str(site.molecule.number),
            site.molecule.name,
            str(site.charge) if site.charge else "0.00",
        ]
        formattedStr = "\t".join(lineList) + "\n"
        out_file.writelines(formattedStr)


def _write_bonds_info(top, out_file):
    """writes the bonds in a topology with assigned atom number."""
    # TODO: need to add bond type info as well
    # TODO: double vs single bonds?
    out_file.write("@<TRIPOS>BOND\n")
    atom_id = {site: idx + 1 for idx, site in enumerate(top.sites)}

    for index, bond in enumerate(top.bonds):
        idx1 = atom_id[bond.connection_members[0]]
        idx2 = atom_id[bond.connection_members[1]]

        idx1, idx2 = min(idx1, idx2), max(idx1, idx2)

        bond_info_str = f"{index + 1} {idx1} {idx2} 1\n"
        out_file.write(bond_info_str)


def _write_molecule_info(top, out_file):
    # TODO: This should check the last three integers to write them correctly
    # TODO: This should handle different molecule types
    out_file.write("@<TRIPOS>MOLECULE\n")
    name = top.name
    num_atoms = top.n_sites
    num_bonds = top.n_bonds
    out_file.write(f"{name}\n{num_atoms} {num_bonds} 1 0 0\n")


def _write_box_info(top, out_file):
    # @<TRIPOS>SUBSTRUCTURE
    # 1 ****        1 TEMP                        0 ****  **** 0 ROOT
    pass
