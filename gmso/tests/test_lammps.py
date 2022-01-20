import pytest
import unyt as u
from unyt.testing import assert_allclose_units

import gmso
from gmso.core.box import Box
from gmso.formats.lammpsdata import read_lammpsdata, write_lammpsdata
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path


class TestLammpsWriter(BaseTest):
    @pytest.mark.parametrize(
        "fname", ["data.lammps", "data.data", "data.lammpsdata"]
    )
    def test_write_lammps(self, fname, typed_ar_system):
        print(fname)
        typed_ar_system.save(fname)

    def test_write_lammps_triclinic(self, typed_ar_system):
        typed_ar_system.box = Box(lengths=[1, 1, 1], angles=[60, 90, 120])
        typed_ar_system.save("triclinic.lammps")

    def test_ethane_lammps(self, typed_ethane):
        typed_ethane.save("ethane.lammps")

    def test_water_lammps(self, typed_water_system):
        typed_water_system.save("data.lammps")

    def test_read_lammps(self, filename=get_path("data.lammps")):
        top = gmso.Topology.load(filename)

    def test_read_box(self, filename=get_path("data.lammps")):
        read = gmso.Topology.load(filename)

        assert read.box == Box(lengths=[1, 1, 1])

    def test_read_n_sites(self, typed_ar_system):
        typed_ar_system.save("ar.lammps")
        read = gmso.Topology.load("ar.lammps")

        assert read.n_sites == 100

    def test_read_mass(self, filename=get_path("data.lammps")):
        read = gmso.Topology.load(filename)
        masses = [i.mass for i in read.atom_types]

        assert_allclose_units(
            masses, u.unyt_array(1.0079, u.g / u.mol), rtol=1e-5, atol=1e-8
        )

    def test_read_charge(self, filename=get_path("data.lammps")):
        read = gmso.Topology.load(filename)
        charge = [i.charge for i in read.atom_types]

        assert_allclose_units(
            charge, u.unyt_array(0, u.C), rtol=1e-5, atol=1e-8
        )

    def test_read_sigma(self, filename=get_path("data.lammps")):
        read = gmso.Topology.load(filename)
        lj = [i.parameters for i in read.atom_types][0]

        assert_allclose_units(
            lj["sigma"], u.unyt_array(3, u.angstrom), rtol=1e-5, atol=1e-8
        )

    def test_read_epsilon(self, filename=get_path("data.lammps")):
        read = gmso.Topology.load(filename)
        lj = [i.parameters for i in read.atom_types][0]

        assert_allclose_units(
            lj["epsilon"],
            u.unyt_array(0.0717, (u.kcal / u.mol)),
            rtol=1e-5,
            atol=1e-8,
        )

    def test_read_water(self, typed_water_system):
        typed_water_system.save("water.lammps")
        water = gmso.Topology.load("water.lammps")

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
        typed_ar_system.save("triclinic.lammps")

        read = gmso.Topology.load("triclinic.lammps")
        assert_allclose_units(
            read.box.lengths,
            u.unyt_array([1, 1, 1], u.nm),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            read.box.angles,
            u.unyt_array([60, 90, 120], u.degree),
            rtol=1e-5,
            atol=1e-8,
        )

    def test_read_n_bonds(self, typed_ethane):
        typed_ethane.save("ethane.lammps")
        read = gmso.Topology.load("ethane.lammps")

        assert read.n_bonds == 7

    def test_read_n_angles(self, typed_ethane):
        typed_ethane.save("ethane.lammps")
        read = gmso.Topology.load("ethane.lammps")

        assert read.n_angles == 12

    def test_read_bond_params(self, typed_ethane):
        typed_ethane.save("ethane.lammps")
        read = gmso.Topology.load("ethane.lammps")
        bond_params = [i.parameters for i in read.bond_types]

        assert_allclose_units(
            bond_params[0]["k"],
            u.unyt_array(680, (u.kcal / u.mol / u.angstrom / u.angstrom)),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            bond_params[0]["r_eq"],
            u.unyt_array(1.09, (u.angstrom)),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            bond_params[1]["k"],
            u.unyt_array(536, (u.kcal / u.mol / u.angstrom / u.angstrom)),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            bond_params[1]["r_eq"],
            u.unyt_array(1.529, (u.angstrom)),
            rtol=1e-5,
            atol=1e-8,
        )

    def test_read_angle_params(self, typed_ethane):
        typed_ethane.save("ethane.lammps")
        read = gmso.Topology.load("ethane.lammps")
        angle_params = [i.parameters for i in read.angle_types]

        assert_allclose_units(
            angle_params[0]["k"],
            u.unyt_array(75, (u.kcal / u.mol / u.radian / u.radian)),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            angle_params[0]["theta_eq"],
            u.unyt_array(110.7, (u.degree)),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            angle_params[1]["k"],
            u.unyt_array(66, (u.kcal / u.mol / u.radian / u.radian)),
            rtol=1e-5,
            atol=1e-8,
        )
        assert_allclose_units(
            angle_params[1]["theta_eq"],
            u.unyt_array(107.8, (u.degree)),
            rtol=1e-5,
            atol=1e-8,
        )


"""
    def test_read_n_diherals(self, typed_ethane):
        typed_ethane.save("ethane.lammps")
        read = gmso.Topology.load("ethane.lammps")

        assert read.n_dihedrals == 9
"""
