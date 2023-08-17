import numpy as np
import pytest
import unyt as u

from gmso.exceptions import EngineIncompatibilityError
from gmso.formats.mcf import write_mcf
from gmso.tests.base_test import BaseTest


def parse_mcf(filename):
    mcf_data = []
    mcf_idx = {}
    with open(filename) as f:
        for line in f:
            mcf_data.append(line.strip().split())

    for idx, line in enumerate(mcf_data):
        if len(line) > 1:
            if line[1] == "Atom_Info":
                mcf_idx["Atom_Info"] = idx
        if len(line) > 1:
            if line[1] == "Bond_Info":
                mcf_idx["Bond_Info"] = idx
        if len(line) > 1:
            if line[1] == "Angle_Info":
                mcf_idx["Angle_Info"] = idx
        if len(line) > 1:
            if line[1] == "Dihedral_Info":
                mcf_idx["Dihedral_Info"] = idx
        if len(line) > 1:
            if line[1] == "Fragment_Info":
                mcf_idx["Fragment_Info"] = idx
        if len(line) > 1:
            if line[1] == "Fragment_Connectivity":
                mcf_idx["Fragment_Connectivity"] = idx

    return mcf_data, mcf_idx


def is_charge_neutral(mcf_data, mcf_idx):
    n_sites = int(mcf_data[mcf_idx["Atom_Info"] + 1][0])
    net_q = 0.0
    for idx, site in enumerate(range(n_sites)):
        net_q += float(mcf_data[mcf_idx["Atom_Info"] + 2 + idx][4])
    return np.isclose(
        net_q,
        0.0,
    )


class TestMCF(BaseTest):
    def test_write_lj_full(self, n_typed_ar_system):
        top = n_typed_ar_system(n_sites=1)
        top.save("ar.mcf")

        mcf_data, mcf_idx = parse_mcf("ar.mcf")

        assert is_charge_neutral(mcf_data, mcf_idx)

        assert mcf_data[mcf_idx["Atom_Info"] + 1][0] == "1"
        assert mcf_data[mcf_idx["Atom_Info"] + 2][1] == "Ar"
        assert mcf_data[mcf_idx["Atom_Info"] + 2][2] == "Ar"
        assert mcf_data[mcf_idx["Atom_Info"] + 2][5] == "LJ"
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][3]),
            top.sites[0].mass.in_units(u.amu).value,
        )
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][4]),
            top.sites[0].charge.in_units(u.elementary_charge).value,
        )
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][6]),
            (top.sites[0].atom_type.parameters["epsilon"] / u.kb)
            .in_units(u.K)
            .value,
        )
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][7]),
            top.sites[0]
            .atom_type.parameters["sigma"]
            .in_units(u.Angstrom)
            .value,
        )
        assert mcf_data[mcf_idx["Bond_Info"] + 1][0] == "0"
        assert mcf_data[mcf_idx["Angle_Info"] + 1][0] == "0"
        assert mcf_data[mcf_idx["Dihedral_Info"] + 1][0] == "0"

        assert mcf_data[mcf_idx["Fragment_Info"] + 1][0] == "1"
        assert mcf_data[mcf_idx["Fragment_Info"] + 2] == ["1", "1", "1"]

        assert mcf_data[mcf_idx["Fragment_Connectivity"] + 1][0] == "0"

        assert np.allclose(float(mcf_data[-5][0]), 0.0)
        assert np.allclose(float(mcf_data[-5][1]), 0.0)
        assert np.allclose(float(mcf_data[-5][2]), 0.5)
        assert np.allclose(float(mcf_data[-5][3]), 1.0)
        assert np.allclose(float(mcf_data[-4][0]), 0.0)
        assert np.allclose(float(mcf_data[-4][1]), 0.0)
        assert np.allclose(float(mcf_data[-4][2]), 0.5)
        assert np.allclose(float(mcf_data[-4][3]), 1.0)

    def test_write_not_neutral(self, n_typed_ar_system):
        top = n_typed_ar_system(n_sites=1)
        top.sites[0].charge = 1.0 * u.elementary_charge
        with pytest.raises(ValueError):
            top.save("ar.mcf")

    def test_write_mie_full(self, n_typed_xe_mie):
        top = n_typed_xe_mie()
        top.save("xe.mcf")

        mcf_data, mcf_idx = parse_mcf("xe.mcf")

        assert is_charge_neutral(mcf_data, mcf_idx)
        # Check a some atom info
        assert mcf_data[mcf_idx["Atom_Info"] + 1][0] == "1"
        assert mcf_data[mcf_idx["Atom_Info"] + 2][1] == "Xe"
        assert mcf_data[mcf_idx["Atom_Info"] + 2][2] == "Xe"
        assert mcf_data[mcf_idx["Atom_Info"] + 2][5] == "Mie"
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][3]),
            top.sites[0].mass.in_units(u.amu).value,
        )
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][4]),
            top.sites[0].charge.in_units(u.elementary_charge).value,
        )
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][6]),
            (top.sites[0].atom_type.parameters["epsilon"] / u.kb)
            .in_units(u.K)
            .value,
        )
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][7]),
            top.sites[0]
            .atom_type.parameters["sigma"]
            .in_units(u.Angstrom)
            .value,
        )
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][8]),
            top.sites[0].atom_type.parameters["n"],
        )
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][9]),
            top.sites[0].atom_type.parameters["m"],
        )

        assert mcf_data[mcf_idx["Bond_Info"] + 1][0] == "0"
        assert mcf_data[mcf_idx["Angle_Info"] + 1][0] == "0"
        assert mcf_data[mcf_idx["Dihedral_Info"] + 1][0] == "0"

        assert mcf_data[mcf_idx["Fragment_Info"] + 1][0] == "1"
        assert mcf_data[mcf_idx["Fragment_Info"] + 2] == ["1", "1", "1"]

        assert mcf_data[mcf_idx["Fragment_Connectivity"] + 1][0] == "0"

        assert np.allclose(float(mcf_data[-5][0]), 0.0)
        assert np.allclose(float(mcf_data[-5][1]), 0.0)
        assert np.allclose(float(mcf_data[-5][2]), 0.5)
        assert np.allclose(float(mcf_data[-5][3]), 1.0)
        assert np.allclose(float(mcf_data[-4][0]), 0.0)
        assert np.allclose(float(mcf_data[-4][1]), 0.0)
        assert np.allclose(float(mcf_data[-4][2]), 0.5)
        assert np.allclose(float(mcf_data[-4][3]), 1.0)

    def test_modified_incompatible_expressions(self, typed_ethane):
        top = typed_ethane

        # Test that we can't write a MCF with a modified potential
        next(iter(top.atom_types)).set_expression("sigma + epsilon*r")
        with pytest.raises(EngineIncompatibilityError):
            top.save("lj.mcf")

        next(iter(top.bond_types)).set_expression("k*(r-r_eq)**3")
        with pytest.raises(EngineIncompatibilityError):
            top.save("bond.mcf")

        # Modified angles
        next(iter(top.angle_types)).set_expression(
            "0.5 * k*(theta-theta_eq)**3"
        )
        with pytest.raises(EngineIncompatibilityError):
            top.save("angle.mcf")

        # Modified dihedrals
        # next(iter(top.dihedral_types)).set_expression("c0 * c1 * c2 * c3 * c4 * c5")
        # with pytest.raises(EngineIncompatibilityError):
        #    top.save("dihedral.mcf")

    def test_typed_ethane(self, typed_ethane):
        top = typed_ethane
        write_mcf(top, "ethane.mcf")

        mcf_data, mcf_idx = parse_mcf("ethane.mcf")

        assert is_charge_neutral(mcf_data, mcf_idx)

        # Check atom info
        assert mcf_data[mcf_idx["Atom_Info"] + 1][0] == "8"
        assert mcf_data[mcf_idx["Atom_Info"] + 2][1] == "opls_135"
        assert mcf_data[mcf_idx["Atom_Info"] + 2][2] == "C"
        assert mcf_data[mcf_idx["Atom_Info"] + 2][5] == "LJ"
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][3]),
            top.sites[0].mass.in_units(u.amu).value,
        )
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][4]),
            top.sites[0].charge.in_units(u.elementary_charge).value,
        )
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][6]),
            (top.sites[0].atom_type.parameters["epsilon"] / u.kb)
            .in_units(u.K)
            .value,
        )
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][7]),
            top.sites[0]
            .atom_type.parameters["sigma"]
            .in_units(u.Angstrom)
            .value,
        )

        # Check bond section
        assert mcf_data[mcf_idx["Bond_Info"] + 1][0] == "7"
        assert mcf_data[mcf_idx["Bond_Info"] + 2][3] == "fixed"
        assert np.isclose(
            float(mcf_data[mcf_idx["Bond_Info"] + 2][4]),
            top.bonds[0]
            .bond_type.parameters["r_eq"]
            .in_units(u.Angstrom)
            .value,
        )

        # Check angle section
        assert mcf_data[mcf_idx["Angle_Info"] + 1][0] == "12"
        assert mcf_data[mcf_idx["Angle_Info"] + 2][4] == "harmonic"
        assert np.isclose(
            float(mcf_data[mcf_idx["Angle_Info"] + 2][6]),
            top.angles[0]
            .angle_type.parameters["theta_eq"]
            .in_units(u.degree)
            .value,
        )
        # TODO: account for 0.5 factors
        # assert np.isclose(
        # float(mcf_data[mcf_idx["Angle_Info"] + 2][5]),
        # (top.angles[0].angle_type.parameters["k"] / u.kb)
        # .in_units(u.K / u.radian**2)
        # .value,
        # )

        # Check dihedral section
        assert mcf_data[mcf_idx["Dihedral_Info"] + 1][0] == "9"
        dihedral_style = mcf_data[mcf_idx["Dihedral_Info"] + 2][5].lower()
        assert dihedral_style.lower() == "opls"

        # TODO: account for 0.5 factors
        # for idx, k in enumerate(["k1", "k2", "k3", "k4"]):
        # assert np.isclose(
        # float(mcf_data[mcf_idx["Dihedral_Info"] + 2][6 + idx]),
        # top.dihedrals[0]
        # .dihedral_type.parameters[k]
        # .in_units(u.kilojoule / u.mole)
        # .value,
        # )
