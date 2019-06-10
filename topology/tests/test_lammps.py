from topology.core.box import Box
from topology.formats.lammpsdata import write_lammpsdata, read_lammpsdata
from topology.tests.base_test import BaseTest


class TestLammpsWriter(BaseTest):
    def test_write_lammps(self, topology_site):
        top = topology_site()
        write_lammpsdata(top, filename='data.lammps')

    def test_write_lammps_triclinic(self, topology_site):
        top = topology_site()
        top.box = Box(lengths=[1,1,1], angles=[60,90,120])
        write_lammpsdata(top, filename='data.triclinic')

    def test_read_lammps(self, topology_site):
        top = topology_site()
        write_lammpsdata(top, filename='data.lammps')
        read_lammpsdata('data.lammps')

    def test_read_sites(self, topology_site):
        top = topology_site()
        write_lammpsdata(top, filename='data.lammps')
        read = read_lammpsdata('data.lammps')

        assert read.box == Box(lengths=[1, 1, 1])

    def test_read_n_sites(self, topology_site):
        top = topology_site(sites=4)
        write_lammpsdata(top, filename='data.lammps')
        read = read_lammpsdata('data.lammps')

        assert read.n_sites == 4
