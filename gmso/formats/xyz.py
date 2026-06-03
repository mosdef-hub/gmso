"""Read and write XYZ files."""

import datetime
from pathlib import Path
from typing import Union

import numpy as np
import unyt as u

from gmso.core.atom import Atom
from gmso.core.topology import Topology
from gmso.formats.formats_registry import loads_as, saves_as


@loads_as(".xyz")
def read_xyz(filename: Union[str, Path]) -> Topology:
    """Read an XYZ file and return a :class:`~gmso.Topology`.

    Parameters
    ----------
    filename : str or pathlib.Path
        Path to the ``.xyz`` file to read.

    Returns
    -------
    gmso.Topology
        Topology containing the sites parsed from *filename*.  No bonds or
        box information are set; call :meth:`~gmso.Topology.update_topology`
        afterwards if needed.

    Raises
    ------
    ValueError
        If the number of coordinate lines does not match the atom count in
        the first line of the file.
    """
    top = Topology()

    with open(filename, "r") as xyz_file:
        n_atoms = int(xyz_file.readline())
        xyz_file.readline()
        coords = np.zeros(shape=(n_atoms, 3)) * u.nanometer
        for row, _ in enumerate(coords):
            line = xyz_file.readline().split()
            if not line:
                msg = (
                    "Incorrect number of lines in input file. Based on the "
                    "number in the first line of the file, {} rows of atoms "
                    "were expected, but at least one fewer was found."
                )
                raise ValueError(msg.format(n_atoms))
            tmp = np.array(line[1:4], dtype=float) * u.angstrom
            coords[row] = tmp.in_units(u.nanometer)
            site = Atom(name=line[0], position=coords[row])
            top.add_site(site)
        top.update_topology()

        line = xyz_file.readline().split()
        if line:
            msg = (
                "Incorrect number of lines in input file. Based on the "
                "number in the first line of the file, {} rows of atoms "
                "were expected, but at least one more was found."
            )
            raise ValueError(msg.format(n_atoms))

    return top


@saves_as(".xyz")
def write_xyz(
    top: Topology,
    filename: Union[str, Path],
    decimals: int = 3,
) -> None:
    """Write a :class:`~gmso.Topology` to an XYZ file.

    Parameters
    ----------
    top : gmso.Topology
        Topology to write.
    filename : str or pathlib.Path
        Destination file path.
    decimals : int, optional, default=3
        Number of decimal places written for each coordinate value.
    """
    with open(filename, "w") as out_file:
        out_file.write("{:d}\n".format(top.n_sites))
        out_file.write(
            "{} {} written by topology at {}\n".format(
                top.name, filename, str(datetime.datetime.now())
            )
        )
        out_file.write(_prepare_particles(top, decimals))


def _prepare_particles(top: Topology, decimals: int) -> str:
    atom_info = str()
    for _, site in enumerate(top.sites):
        if site.element is not None:
            tmp_name = site.element.symbol
        else:
            tmp_name = "X"

        x = site.position[0].in_units(u.angstrom).value
        y = site.position[1].in_units(u.angstrom).value
        z = site.position[2].in_units(u.angstrom).value
        atom_info = (
            atom_info
            + f"{tmp_name} {x:{decimals + 5}.{decimals}f} {y:{decimals + 5}.{decimals}f} {z:{decimals + 5}.{decimals}f}\n"
        )
    return atom_info
