import unyt as u
from unyt.testing import assert_allclose_units

import gmso
from gmso.core.box import Box
from gmso.formats.lammpsdata import read_lammpsdata, write_lammpsdata
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path


class TestLammpsWriter(BaseTest):
    def test_write_lammps(self, typed_ar_system):
        write_lammpsdata(typed_ar_system, filename="data.lammps")

    def test_write_lammps_triclinic(self, typed_ar_system):
        typed_ar_system.box = Box(lengths=[1, 1, 1], angles=[60, 90, 120])
        write_lammpsdata(typed_ar_system, filename="data.triclinic")

    def test_ethane_lammps(self, typed_ethane):
        write_lammpsdata(typed_ethane, "data.ethane")

    def test_water_lammps(self, typed_water_system):
        write_lammpsdata(typed_water_system, "data.water")

    def test_read_lammps(self, filename=get_path("data.lammps")):
        read_lammpsdata(filename)

    def test_read_sites(self, filename=get_path("data.lammps")):
        read = read_lammpsdata(filename)

        assert read.box == Box(lengths=[1, 1, 1])

    def test_read_n_sites(self, typed_ar_system):
        write_lammpsdata(typed_ar_system, filename="data.ar")
        read = read_lammpsdata("data.ar")

        assert read.n_sites == 100

    def test_read_mass(self, filename=get_path("data.lammps")):
        read = read_lammpsdata(filename)
        masses = [i.mass for i in read.atom_types]

        assert_allclose_units(
            masses, u.unyt_array(1.0079, u.g), rtol=1e-5, atol=1e-8
        )

    def test_read_charge(self, filename=get_path("data.lammps")):
        read = read_lammpsdata(filename)
        charge = [i.charge for i in read.atom_types]

        assert_allclose_units(
            charge, u.unyt_array(0, u.C), rtol=1e-5, atol=1e-8
        )

    def test_read_sigma(self, filename=get_path("data.lammps")):
        read = read_lammpsdata(filename)
        lj = [i.parameters for i in read.atom_types][0]

        assert_allclose_units(
            lj["sigma"], u.unyt_array(3, u.angstrom), rtol=1e-5, atol=1e-8
        )

    def test_read_epsilon(self, filename=get_path("data.lammps")):
        read = read_lammpsdata(filename)
        lj = [i.parameters for i in read.atom_types][0]

        assert_allclose_units(
            lj["epsilon"],
            u.unyt_array(0.0717, (u.kcal / u.mol)),
            rtol=1e-5,
            atol=1e-8,
        )

    def test_read_water(self, typed_water_system):
        write_lammpsdata(typed_water_system, filename="data.water")
        water = read_lammpsdata("data.water")

        assert_allclose_units(
            water.sites[0].charge,
            u.unyt_array(-0.834, u.elementary_charge),
            rtol=1e-5,
            atol=1e-8,
        )
        assert water.n_sites == 6
        assert water.n_connections == 6

    def test_read_lammps_triclinic(self, typed_ar_system):
        typed_ar_system.box = Box(lengths=[1, 1, 1], angles=[60, 90, 120])
        write_lammpsdata(typed_ar_system, filename="data.triclinic")

        read = read_lammpsdata("data.triclinic")
        assert u.allclose_units(
            read.box.lengths, u.unyt_array([1.0, 1.0, 1.0], u.nm)
        )
        assert u.allclose_units(
            read.box.angles, u.unyt_array([60, 90, 120], u.degree)
        )
