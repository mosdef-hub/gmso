"""Read and write Gromos87 (.GRO) file format."""

import datetime
import logging
from pathlib import Path
from typing import Union

import numpy as np
import unyt as u
from unyt.array import allclose_units

import gmso
from gmso.core.atom import Atom
from gmso.core.box import Box
from gmso.core.topology import Topology
from gmso.formats.formats_registry import loads_as, saves_as

logger = logging.getLogger(__name__)


@loads_as(".gro")
def read_gro(filename: Union[str, Path]) -> Topology:
    """Read a Gromos87 (GRO) file and return a :class:`~gmso.Topology`.

    The GRO format is a plain-text structure file used by the GROMACS
    simulation engine.  It encodes simulation-box parameters, atom count,
    residue/atom names, positions, and (optionally) velocities.

    Parameters
    ----------
    filename : str or pathlib.Path
        Path to the ``.gro`` file to read.

    Returns
    -------
    gmso.Topology
        Topology containing the sites and box parsed from *filename*.

    Notes
    -----
    * GRO files do not specify bonds; the returned topology will have no
      connections.
    * Velocities are detected but not parsed.
    * Only single-frame GRO files are supported.
    * Residue/resid information is mapped to ``site.molecule`` and
      ``site.residue`` but the original resid integers are rebased to
      0-indexed integers.
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
            res_id = int(line[:5].strip()) - 1  # reformat from 1 to 0 index in gmso
            res_name = line[5:10].strip()
            atom_name = line[10:15].strip()

            positions = line[20:].split()
            coords[row] = u.nm * np.array(
                [
                    float(positions[0]),
                    float(positions[1]),
                    float(positions[2]),
                ]
            )
            site = Atom(name=atom_name, position=coords[row])

            site.molecule = (res_name, res_id)
            site.residue = (res_name, res_id)
            top.add_site(site, update_types=False)

        if len(positions) == 6:
            logger.info("Velocity information presents but will not be parsed.")
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
def write_gro(
    top: Topology,
    filename: Union[str, Path],
    n_decimals: int = 3,
    shift_coord: bool = False,
) -> None:
    """Write a :class:`~gmso.Topology` to a Gromos87 (GRO) file.

    Parameters
    ----------
    top : gmso.Topology
        Topology to write.
    filename : str or pathlib.Path
        Destination file path.
    n_decimals : int, optional, default=3
        Number of decimal places to use for coordinate values.
    shift_coord : bool, optional, default=False
        When ``True``, shift all coordinates so that the minimum position in
        each dimension is zero.  Useful for visualisation; not required by the
        GRO format.

    Notes
    -----
    * Velocities are not written.
    * Residue/resid indices cycle back at 99 999 to comply with the fixed-width
      GRO column format.
    """
    pos_array = np.ndarray.copy(top.positions)
    if shift_coord:
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
        out_file.write(_prepare_atoms(top, pos_array, n_decimals))
        out_file.write(_prepare_box(top))


def _validate_positions(pos_array):
    """Modify coordinates, as necessary, to fit limitations of the GRO format."""
    if np.min(pos_array) < 0:
        logger.info(
            "Topology contains some negative positions. Translating "
            "in order to ensure all coordinates are non-negative."
        )
    min_xyz = np.min(pos_array, axis=0)
    min_xyz0 = np.where(min_xyz < 0 * min_xyz.units, min_xyz, 0 * min_xyz.units)

    pos_array -= min_xyz0

    return pos_array


def _prepare_atoms(top, updated_positions, n_decimals):
    out_str = str()
    logger.info(
        "Residue information is parsed from site.molecule,"
        "or site.residue if site.molecule does not exist."
        "Note that the residue idx will be bumped by 1 since GROMACS utilize 1-index."
    )
    # we need to sort through the sites to provide a unique number for each molecule/residue
    # we will store the unique id in dictionary where the key is the idx
    site_res_id = dict()
    seen = dict()
    for idx, site in enumerate(top.sites):
        if site.molecule:
            if site.molecule not in seen:
                seen[site.molecule] = len(seen) + 1
            site_res_id[idx] = seen[site.molecule]
        elif site.residue:
            if site.residue not in seen:
                seen[site.residue] = len(seen) + 1
            site_res_id[idx] = seen[site.residue]
        else:
            if "MOL" not in seen:
                seen["MOL"] = len(seen) + 1
            site_res_id[idx] = seen["MOL"]

    for idx, (site, pos) in enumerate(zip(top.sites, updated_positions)):
        if site.molecule:
            res_id = site_res_id[idx]
            res_name = (
                site.molecule.name
                if len(site.molecule.name) <= 5
                else site.molecule.name[:5]
            )

            site.label = f"res_id: {res_id}, " + site.label
        elif site.residue:
            res_id = site_res_id[idx]
            res_name = (
                site.residue.name
                if len(site.residue.name) <= 5
                else site.residue.name[:5]
            )
            site.label = f"res_id: {res_id}, " + site.label
        else:
            res_id = site_res_id[idx]
            res_name = "MOL"

            site.label = f"res_id: {res_id}, " + site.label

        atom_name = site.name if len(site.name) <= 5 else site.name[:5]
        atom_id = idx + 1

        # gromacs doesn't actually use the atom id in the .gro file
        # so we will just loop back to 1 once we exceed 99999
        # as is suggested in the FAQ in the manual.

        max_val = 99999
        atom_id = atom_id % max_val
        res_id = res_id % max_val

        varwidth = 5 + n_decimals
        crdfmt = f"{{:{varwidth}.{n_decimals}f}}"

        # preformat pos str
        crt_x = crdfmt.format(pos[0].in_units(u.nm).value)[:varwidth]
        crt_y = crdfmt.format(pos[1].in_units(u.nm).value)[:varwidth]
        crt_z = crdfmt.format(pos[2].in_units(u.nm).value)[:varwidth]
        out_str = out_str + "{0:5d}{1:5s}{2:>5s}{3:5d}{4}{5}{6}\n".format(
            res_id,
            res_name,
            atom_name,
            atom_id,
            crt_x,
            crt_y,
            crt_z,
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
        out_str = out_str + " {:0.5f} {:0.5f} {:0.5f}\n".format(
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
