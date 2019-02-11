from __future__ import division

from collections import OrderedDict
from warnings import warn
import itertools as it

import numpy as np


from topology.core.topology import Topology
from topology.core.site import Site


def write_lammpsdata(topology, filename, atom_style='full',
        nbfix_in_data_file=True):
    """Output a LAMMPS data file.
    
    Outputs a LAMMPS data file in the 'full' atom style format. Assumes use
    of 'real' units. See http://lammps.sandia.gov/doc/atom_style.html for
    more information on atom styles.

    Parameters
    ----------
    Topology : Topology
        MoSDeF Topology Object
    filename : str
        Path of the output file
    atom_style: str
        Defines the style of atoms to be saved in a LAMMPS data file. The following atom
        styles are currently supported: 'full', 'atomic', 'charge', 'molecular'
        see http://lammps.sandia.gov/doc/atom_style.html for more
        information on atom styles.

    Notes
    -----
    See http://lammps.sandia.gov/doc/2001/data_format.html for a full description
    of the LAMMPS data format. Currently the following sections are supported (in
    addition to the header): *Masses*, *Nonbond Coeffs*, *Bond Coeffs*, *Angle
    Coeffs*, *Dihedral Coeffs*, *Atoms*, *Bonds*, *Angles*, *Dihedrals*

    Some of this function has beed adopted from `mdtraj`'s support of the LAMMPSTRJ
    trajectory format. See https://github.com/mdtraj/mdtraj/blob/master/mdtraj/formats/lammpstrj.py for details.

    """
    if atom_style not in ['atomic', 'charge', 'molecular', 'full']:
        raise ValueError('Atom style "{}" is invalid or is not currently supported'.format(atom_style))

    xyz = np.array([[site.xx,site.xy,site.xz] for site in topology.sites])

    forcefield = True
    if structure[0].type == '':
        forcefield = False

    # Internally use nm
    box = Box(lengths=np.array([0.1 * val for val in topology.box.lengths]),
              angles=topology.box.angles)

    if forcefield:
        types = [site.type for site in topology.sites]
    else:
        types = [site.name for site in topology.sites]

    unique_types = list(set(types))
    unique_types.sort(key=natural_sort)

    charges = [site.charge for site in structure.sites]
