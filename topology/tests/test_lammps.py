from topology.core.box import Box
from topology.formats.lammpsdata import write_lammpsdata, read_lammpsdata
from topology.tests.base_test import BaseTest
import unyt as u


class TestLammpsWriter(BaseTest):
    def test_write_lammps(self, topology_site):
        top = topology_site()
        write_lammpsdata(top, filename='data.lammps')

    def test_write_lammps_triclinic(self, topology_site):
        topology_site().box = Box(lengths=[1,1,1], angles=[60,90,120])
        write_lammpsdata(topology_site(), filename='data.triclinic')

    def test_read_lammps(self, topology_site):
        write_lammpsdata(topology_site(), filename='data.lammps')
        read_lammpsdata('data.lammps')

    def test_read_sites(self, topology_site):
        write_lammpsdata(topology_site(), filename='data.lammps')
        read = read_lammpsdata('data.lammps')

        assert read.box == Box(lengths=[1, 1, 1])

    def test_read_n_sites(self, topology_site):
        write_lammpsdata(topology_site(sites=4),
                filename='data.lammps')
        read = read_lammpsdata('data.lammps')

        assert read.n_sites == 4

    def test_read_mass(self, topology_site):
        write_lammpsdata(topology_site(), filename='data.lammps')
        read = read_lammpsdata('data.lammps')
        masses = [i.mass for i in read.atom_types]

        assert masses == u.unyt_array(1, u.g)

    def test_read_mass(self, topology_site):
        write_lammpsdata(topology_site(), filename='data.lammps')
        read = read_lammpsdata('data.lammps')
        charge = [i.charge for i in read.atom_types]

        assert charge == u.unyt_array(0, u.C)

    def test_read_sigma(self, topology_site):
        write_lammpsdata(topology_site(), filename='data.lammps')
        read = read_lammpsdata('data.lammps')
        lj = [i.parameters for i in read.atom_types][0]

        assert lj['sigma'] == u.unyt_array(3, u.angstrom)

    def test_read_epsilon(self, topology_site):
        write_lammpsdata(topology_site(), filename='data.lammps')
        read = read_lammpsdata('data.lammps')
        lj = [i.parameters for i in read.atom_types][0]

        assert lj['epsilon'] == u.unyt_array(0.0717, (u.kcal/u.mol))
