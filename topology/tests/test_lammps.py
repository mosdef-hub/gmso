import numpy as np
import pytest 
import unyt as u

from topology.core.topology import Topology
from topology.core.site import Site
from topology.core.bond import Bond
from topology.core.box import Box
from topology.testing.utils import allclose
from topology.formats.lammpsdata import write_lammpsdata
from topology.core.atom_type import AtomType
#import MDAnalysis as mda
#from MDAnalysis import Universe


def test_write_lammps():
    top = Topology()
    top.box = Box(lengths=[1,1,1])
    site1 = Site(name='site1', atom_type=AtomType(mass=1))
    site2 = Site(name='site2', atom_type=AtomType(mass=1))
    connect = Bond(connection_members=[site1, site2])

    top.add_site(site1)
    top.add_site(site2)
    top.add_connection(connect)

    write_lammpsdata(top, filename='data.lammps')

def test_write_lammps_triclinic():
    top = Topology()
    top.box = Box(lengths=[1,1,1], angles=[60,90,120])
    site1 = Site(name='site1', atom_type=AtomType(mass=1))
    site2 = Site(name='site2', atom_type=AtomType(mass=1))
    connect = Bond(connection_members=[site1, site2])

    top.add_site(site1)
    top.add_site(site2)
    top.add_connection(connect)

    write_lammpsdata(top, filename='data.triclinic')

#def test_num_atoms():
#    u = mda.Universe('lammps.data')
#
#    assert len(u.atoms) == 2
