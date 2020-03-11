from gmso.core.box import Box
from gmso.formats.lammpsdata import write_lammpsdata, read_lammpsdata
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path
import unyt as u
import numpy as np


class TestLammpsWriter(BaseTest):
    def test_write_lammps(self, topology_site):
        write_lammpsdata(topology_site(), filename='data.lammps')

    def test_write_lammps_triclinic(self, topology_site):
        topology_site().box = Box(lengths=[1,1,1], angles=[60,90,120])
        write_lammpsdata(topology_site(), filename='data.triclinic')

    def test_read_lammps(self, filename=get_path('data.lammps')):
        read_lammpsdata(filename)

    def test_read_sites(self, filename=get_path('data.lammps')):
        read = read_lammpsdata(filename)

        assert read.box == Box(lengths=[1, 1, 1])

    def test_read_n_sites(self, topology_site):
        write_lammpsdata(topology_site(sites=4),
                filename='data.four_sites')
        read = read_lammpsdata('data.four_sites')

        assert read.n_sites == 4

    def test_read_mass(self, filename=get_path('data.lammps')):
        read = read_lammpsdata(filename)
        masses = [i.mass for i in read.atom_types]

        assert masses == u.unyt_array(1, u.g)

    def test_read_mass(self, filename=get_path('data.lammps')):
        read = read_lammpsdata(filename)
        charge = [i.charge for i in read.atom_types]

        assert charge == u.unyt_array(0, u.C)

    def test_read_sigma(self, filename=get_path('data.lammps')):
        read = read_lammpsdata(filename)
        lj = [i.parameters for i in read.atom_types][0]

        assert np.allclose(lj['sigma'].value,
                u.unyt_array(3, u.angstrom).value)

    def test_read_epsilon(self, filename=get_path('data.lammps')):
        read = read_lammpsdata(filename)
        lj = [i.parameters for i in read.atom_types][0]

        assert np.allclose(lj['epsilon'].value,
                u.unyt_array(0.0717, (u.kcal/u.mol)).value)
