from gmso.core.box import Box
from gmso.formats.lammpsdata import write_lammpsdata
from gmso.tests.base_test import BaseTest
import unyt as u


class TestLammpsWriter(BaseTest):
    def test_write_lammps(self, topology_site):
        top = topology_site()
        write_lammpsdata(top, filename='data.lammps')

    def test_write_lammps_triclinic(self, topology_site):
        top = topology_site()
        top.box = Box(lengths=[1,1,1], angles=[60,90,120])
        write_lammpsdata(top, filename='data.triclinic')

    def test_write_lammps_charges(self, topology_site):
        top = topology_site()
        for site in top.sites:
            site.charge = 1 * u.elementary_charge

        write_lammpsdata(top, filename='data.lammps')
