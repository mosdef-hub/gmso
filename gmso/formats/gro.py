"""Read and write Gromos87 (.GRO) file format."""
import datetime
import warnings

import numpy as np
import unyt as u
from unyt.array import allclose_units

import gmso
from gmso.core.atom import Atom
from gmso.core.box import Box
from gmso.core.topology import Topology
from gmso.exceptions import NotYetImplementedWarning
from gmso.formats.formats_registry import loads_as, saves_as


@loads_as(".gro")
def read_gro(filename):
    """Create a topology from a provided gro file.

    The Gromos87 (gro) format is a common plain text structure file used
    commonly with the GROMACS simulation engine.  This file contains the
    simulation box parameters, number of atoms, the residue and atom number for
    each atom, as well as their positions and velocities (velocity is
    optional).  This method will receive a filepath representation either as a
    string, or a file object and return a `topology`.

    Parameters
    ----------
    filename : str or file object
        The path to the gro file either as a string, or a file object that
        points to the gro file.

    Returns
    -------
    gmso.core.topology
        A `topology` object containing site information


    Notes
    -----
    Gro files do not specify connections between atoms, the returned topology
    will not have connections between sites either.

    Currently this implementation does not support a gro file with more than 1
    frame.

    All residues and resid information from the gro file are currently lost
    when converting to `topology`.

    """
    top = Topology()

    with open(filename, "r") as gro_file:
        top.name = str(gro_file.readline().strip())
        n_atoms = int(gro_file.readline())
        coords = u.nm * np.zeros(shape=(n_atoms, 3))
        for row, _ in enumerate(coords):
            line = gro_file.readline()
            if not line:
                msg = (
                    "Incorrect number of lines in .gro file. Based on the "
                    "number in the second line of the file, {} rows of"
                    "atoms were expected, but at least one fewer was found."
                )
                raise ValueError(msg.format(n_atoms))
            resid = int(line[:5])
            res_name = line[5:10]
            atom_name = line[10:15]
            atom_id = int(line[15:20])
            coords[row] = u.nm * np.array(
                [
                    float(line[20:28]),
                    float(line[28:36]),
                    float(line[36:44]),
                ]
            )
            site = Atom(name=atom_name, position=coords[row])
            top.add_site(site, update_types=False)
        top.update_topology()

        # Box information
        line = gro_file.readline().split()
        top.box = Box(u.nm * np.array([float(val) for val in line[:3]]))

        # Verify we have read the last line by ensuring the next line in blank
        line = gro_file.readline()
        if line:
            msg = (
                "Incorrect number of lines in input file. Based on the "
                "number in the second line of the file, {} rows of atoms "
                "were expected, but at least one more was found."
            )
            raise ValueError(msg.format(n_atoms))

    return top


@saves_as(".gro")
def write_gro(top, filename):
    """Write a topology to a gro file.

    The Gromos87 (gro) format is a common plain text structure file used
    commonly with the GROMACS simulation engine.  This file contains the
    simulation box parameters, number of atoms, the residue and atom number for
    each atom, as well as their positions and velocities (velocity is
    optional).  This method takes a topology object and a filepath string or
    file object and saves a Gromos87 (gro) to disk.

    Parameters
    ----------
    top : gmso.core.topology
        The `topology` to write out to the gro file.
    filename : str or file object
        The location and name of file to save to disk.


    Notes
    -----
    Multiple residue assignment has not been added, each `site` will belong to
    the same resid of 1 currently.

    Velocities are not written out.

    """
    pos_array = np.ndarray.copy(top.positions)
    pos_array = _validate_positions(pos_array)

    with open(filename, "w") as out_file:
        out_file.write(
            "{} written by GMSO {} at {}\n".format(
                top.name if top.name is not None else "",
                gmso.__version__,
                str(datetime.datetime.now()),
            )
        )
        out_file.write("{:d}\n".format(top.n_sites))
        out_file.write(_prepare_atoms(top, pos_array))
        out_file.write(_prepare_box(top))


def _validate_positions(pos_array):
    """Modify coordinates, as necessary, to fit limitations of the GRO format."""
    if np.min(pos_array) < 0:
        warnings.warn(
            "Topology contains some negative positions. Translating "
            "in order to ensure all coordinates are non-negative."
        )
    min_xyz = np.min(pos_array, axis=0)
    for i, minimum in enumerate(min_xyz):
        if minimum < 0.0:
            for loc in pos_array:
                loc[i] = loc[i] - minimum
    return pos_array


def _prepare_atoms(top, updated_positions):
    out_str = str()
    for idx, (site, pos) in enumerate(zip(top.sites, updated_positions)):
        warnings.warn(
            "Residue information is not currently "
            "stored or written to GRO files.",
            NotYetImplementedWarning,
        )
        # TODO: assign residues
        res_id = 1
        res_name = "X"
        atom_name = site.name
        atom_id = idx + 1
        out_str = (
            out_str
            + "{0:5d}{1:5s}{2:5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}\n".format(
                res_id,
                res_name,
                atom_name,
                atom_id,
                pos[0].in_units(u.nm).value,
                pos[1].in_units(u.nm).value,
                pos[2].in_units(u.nm).value,
            )
        )
    return out_str


def _prepare_box(top):
    out_str = str()
    if allclose_units(
        top.box.angles,
        u.degree * [90, 90, 90],
        rtol=1e-5,
        atol=0.1 * u.degree,
    ):
        out_str = out_str + " {:0.5f} {:0.5f} {:0.5f} \n".format(
            top.box.lengths[0].in_units(u.nm).value.round(6),
            top.box.lengths[1].in_units(u.nm).value.round(6),
            top.box.lengths[2].in_units(u.nm).value.round(6),
        )
    else:
        # TODO: Work around GROMACS's triclinic limitations #30
        vectors = top.box.get_vectors()
        out_str = (
            out_str
            + " {:0.5f} {:0.5f} {:0.5f} {:0.5f} {:0.5f} {:0.5f} {:0.5f} {:0.5f} {:0.5f} \n".format(
                vectors[0, 0].in_units(u.nm).value.round(6),
                vectors[1, 1].in_units(u.nm).value.round(6),
                vectors[2, 2].in_units(u.nm).value.round(6),
                vectors[0, 1].in_units(u.nm).value.round(6),
                vectors[0, 2].in_units(u.nm).value.round(6),
                vectors[1, 0].in_units(u.nm).value.round(6),
                vectors[1, 2].in_units(u.nm).value.round(6),
                vectors[2, 0].in_units(u.nm).value.round(6),
                vectors[2, 1].in_units(u.nm).value.round(6),
            )
        )
    return out_str
