import warnings

import numpy as np
import unyt as u

from topology.core.topology import Topology
from topology.core.site import Site
from topology.core.box import Box


def read_gro(filename):
    top = Topology()

    with open(filename, 'r') as gro_file:
        top.name = str(gro_file.readline().strip())
        n_atoms = int(gro_file.readline())
        coords = u.nm * np.zeros(shape=(n_atoms, 3), dtype=np.float64)
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
            coords[row] = u.nm * np.array([
                float(line[20:28]),
                float(line[28:36]),
                float(line[36:44]),
            ])
            site = Site(name=atom_name, position=coords[row])
            top.add_site(site)

        # Box information
        line = gro_file.readline().split()
        top.box = Box(u.nm * np.array([float(val) for val in line[:3]]))

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

def write_gro(top, filename):
    """Write a topology to a gro file."""

    top = _prepare_topology(top)

    with open(filename, 'w') as out_file:
        import pdb; pdb.set_trace()
        out_file.write('{} written by topology\n'.format(
            top.name if top.name is not None else ''))
        out_file.write('{:d}\n'.format(top.n_sites))
        for idx, site in enumerate(top.site_list):
            res_id = 1
            res_name = 'X'
            atom_name = site.name
            atom_id = idx + 1
            out_file.write('{0:5d}{1:5s}{2:5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}\n'.format(
                res_id,
                res_name,
                atom_name,
                atom_id,
                site.position[0].in_units(u.nm).value,
                site.position[1].in_units(u.nm).value,
                site.position[2].in_units(u.nm).value,
            ))

        out_file.write(' {:0.5f} {:0.5f} {:0.5f} \n'.format(
            top.box.lengths[0].in_units(u.nm).value.round(6),
            top.box.lengths[1].in_units(u.nm).value.round(6),
            top.box.lengths[2].in_units(u.nm).value.round(6),
        ))

def _prepare_topology(top):
    """Modify topology, as necessary, to fit limitations of the GRO format."""
    if np.min(top.positions()) < 0:
        warnings.warn('Topology contains some negative positions. Translating '
                      'in order to ensure all coordinates are non-negative.')
        #top.positions() -= np.min(top.positions(), axis=0)

    return top
