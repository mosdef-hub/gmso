import datetime

import numpy as np
import unyt as u

from gmso.core.topology import Topology
from gmso.core.atom import Atom


def read_xyz(filename):
    """Reader for xyz file format.

    Read in an xyz file at the given path and return a Topology object.

    Parameter
    ----------
    filename : str
        Path to .xyz file that need to be read.

    Return
    ------
    top : topology.Topology
        Topology object
    """
    top = Topology()

    with open(filename, 'r') as xyz_file:
        n_atoms = int(xyz_file.readline())
        xyz_file.readline()
        coords = np.zeros(shape=(n_atoms, 3)) * u.nanometer
        for row, _ in enumerate(coords):
            line = xyz_file.readline().split()
            if not line:
                msg = ('Incorrect number of lines in input file. Based on the '
                       'number in the first line of the file, {} rows of atoms '
                       'were expected, but at least one fewer was found.')
                raise ValueError(msg.format(n_atoms))
            tmp = np.array(line[1:4], dtype=np.float) * u.angstrom
            coords[row] = tmp.in_units(u.nanometer)
            site = Atom(name=line[0], position=coords[row])
            top.add_site(site)
        top.update_topology()

        # Verify we have read the last line by ensuring the next line in blank
        line = xyz_file.readline().split()
        if line:
            msg = ('Incorrect number of lines in input file. Based on the '
                   'number in the first line of the file, {} rows of atoms '
                   'were expected, but at least one more was found.')
            raise ValueError(msg.format(n_atoms))

    return top

def write_xyz(top, filename):
    """Writer for xyz file format.

    Write a Topology object to an xyz file at the given path.

    Paramters
    ---------
    top : topology.Topology
        Topology object that needs to be written out.
    filename : str
        Path to file location.
    """
    with open(filename, 'w') as out_file:
        out_file.write('{:d}\n'.format(top.n_sites))
        out_file.write('{} {} written by topology at {}\n'.format(
            top.name,
            filename,
            str(datetime.datetime.now())))
        for idx, site in enumerate(top.sites):
            # TODO: Better handling of element guessing and site naming
            if site.element is not None:
                tmp_name = site.element.symbol
            else:
                tmp_name = 'X'
            out_file.write('{0} {1:8.3f} {2:8.3f} {3:8.3f}\n'.format(
                tmp_name,
                site.position[0].in_units(u.angstrom).value,
                site.position[1].in_units(u.angstrom).value,
                site.position[2].in_units(u.angstrom).value))
