import numpy as np

from topology.core.topology import Topology
from topology.core.site import Site
from topology.core.box import Box


def read_gro(filename):
    top = Topology()

    with open(filename, 'r') as gro_file:
        top.name = str(gro_file.readline().strip())
        n_atoms = int(gro_file.readline())
        coords = np.zeros(shape=(n_atoms, 3), dtype=np.float64)
        for row, _ in enumerate(coords):
            line = gro_file.readline()
            if not line:
                msg = (
                    'Incorrect number of lines in .gro file. Based on the '
                    'number in the second line of the file, {} rows of'
                    'atoms were expected, but at least one fewer was found.'
                )
                raise ValueError(msg.format(n_atoms))
            resid = int(line[:5])
            res_name = line[5:10]
            atom_name = line[10:15]
            atom_id = int(line[15:20])
            coords[row] = np.array([
                float(line[20:28]),
                float(line[28:36]),
                float(line[36:44]),
            ])
            site = Site(name=atom_name, position=coords[row])
            top.add_site(site)

        # Box information
        line = gro_file.readline().split()
        top.box = Box(np.array([float(val) for val in line[:3]]))

        # Verify we have read the last line by ensuring the next line in blank
        line = gro_file.readline()
        if line:
            msg = (
                'Incorrect number of lines in input file. Based on the '
                'number in the second line of the file, {} rows of atoms '
                'were expected, but at least one more was found.'
            )
            raise ValueError(msg.format(n_atoms))

    return top
