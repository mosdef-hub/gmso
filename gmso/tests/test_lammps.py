import gmso
import unyt as u
from gmso.core.box import Box
from gmso.formats.lammpsdata import write_lammpsdata
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path


class TestLammpsWriter(BaseTest):
    def test_write_lammps(self, topology_site):
        top = topology_site()
        write_lammpsdata(top, filename='data.lammps')

    def test_write_lammps_triclinic(self, topology_site):
        top = topology_site()
        top.box = Box(lengths=[1,1,1], angles=[60,90,120])
        write_lammpsdata(top, filename='data.triclinic')

    def test_water_lammps(self, typed_water_system):
        write_lammpsdata(typed_water_system, 'data.water')
