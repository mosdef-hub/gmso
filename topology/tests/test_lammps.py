import numpy as np
import pytest 
import unyt as u

from topology.core.topology import Topology
from topology.core.site import Site
from topology.core.connection import Connection
from topology.core.box import Box
from topology.testing.utils import allclose
from topology.formats.lammpsdata import write_lammpsdata
from topology.core.atom_type import AtomType
from topology.external.convert_openmm import to_openmm


class TestLammpsWriter(BaseTest)
    def test_write_lammps(self, topology_site):
        top = topology_site()
        write_lammpsdata(top, filename='data.lammps')
    
    def test_write_lammps_triclinic(self, topology_site):
        top = topology_site()
        top.box = Box(lengths=[1,1,1], angles=[60,90,120])
        write_lammpsdata(top, filename='data.triclinic')
