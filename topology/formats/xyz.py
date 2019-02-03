import numpy as np
import unyt as u

from topology.core.topology import Topology
from topology.core.site import Site


def read_xyz(filename):
    top = Topology()

    with open(filename, 'r') as xyz_file:
        n_atoms = int(xyz_file.readline())
        xyz_file.readline()
        coords = np.zeros(shape=(n_atoms, 3), dtype=np.float64) * u.nanometer
        for row, _ in enumerate(coords):
            line = xyz_file.readline().split()
            if not line:
                msg = ('Incorrect number of lines in input file. Based on the '
                       'number in the first line of the file, {} rows of atoms '
                       'were expected, but at least one fewer was found.')
                raise ValueError(msg.format(n_atoms))
            tmp = line[1:4] * u.angstrom
            coords[row] = tmp.convert_to_units(u.nanometer)
            site = Site(name=line[0], position=coords[row])
            top.add_site(site)

        # Verify we have read the last line by ensuring the next line in blank
        line = xyz_file.readline().split()
        if line:
            msg = ('Incorrect number of lines in input file. Based on the '
                   'number in the first line of the file, {} rows of atoms '
                   'were expected, but at least one more was found.')
            raise ValueError(msg.format(n_atoms))

    return top
