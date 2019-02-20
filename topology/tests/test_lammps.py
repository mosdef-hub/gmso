import numpy as np
import pytest
import unyt as u

from topology.core.topology import Topology
from topology.core.site import Site
from topology.core.connection import Connection
from topology.core.box import Box
from topology.testing.utils import allclose
from topology.formats.lammpsdata import write_lammpsdata


#def test_write():
#    top = Topology()
#    top.box = Box(lengths=[1,1,1])
#    site1 = Site(name='site1')
#    site2 = Site(name='site2')
#    connect = Connection(site1=site1, site2=site2)
#
#    top.add_site(site1)
#    top.add_site(site2)
#
#    top.update_connection_list()
#    write_lammpsdata(top, filename='test.lammps')

def test_save():
    pass

def test_save_forcefield():
    pass

def test_save_box():
    pass

def test_save_triclinic_box():
    pass

def test_nbfix():
    pass

def test_save_triclinic_box():
    pass
