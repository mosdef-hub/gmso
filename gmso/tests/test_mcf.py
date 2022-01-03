import numpy as np
import pytest
import unyt as u

from gmso.exceptions import EngineIncompatibilityError
from gmso.formats.mcf import write_mcf
from gmso.tests.base_test import BaseTest


class TestMCF(BaseTest):
    def test_write_lj_simple(self, n_typed_ar_system):
        top = n_typed_ar_system(n_sites=1)
        top.save("ar.mcf")

    def test_write_mie_simple(self, n_typed_xe_mie):
        top = n_typed_xe_mie()
        top.save("xe.mcf")

    def test_write_lj_full(self, n_typed_ar_system):
        top = n_typed_ar_system(n_sites=1)
        top.save("ar.mcf")

        mcf_data = []
        with open("ar.mcf") as f:
            for line in f:
                mcf_data.append(line.strip().split())

        for idx, line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == "Atom_Info":
                    atom_section_start = idx

        assert mcf_data[atom_section_start + 1][0] == "1"
        assert mcf_data[atom_section_start + 2][1] == "Ar"
        assert mcf_data[atom_section_start + 2][2] == "Ar"
        assert mcf_data[atom_section_start + 2][5] == "LJ"
        assert np.isclose(
            float(mcf_data[atom_section_start + 2][3]),
            top.sites[0].mass.in_units(u.amu).value,
        )
        assert np.isclose(
            float(mcf_data[atom_section_start + 2][4]),
            top.sites[0].charge.in_units(u.elementary_charge).value,
        )
        assert np.isclose(
            float(mcf_data[atom_section_start + 2][6]),
            (top.sites[0].atom_type.parameters["epsilon"] / u.kb)
            .in_units(u.K)
            .value,
        )
        assert np.isclose(
            float(mcf_data[atom_section_start + 2][7]),
            top.sites[0]
            .atom_type.parameters["sigma"]
            .in_units(u.Angstrom)
            .value,
        )

    def test_write_mie_full(self, n_typed_xe_mie):
        top = n_typed_xe_mie()
        top.save("xe.mcf")

        mcf_data = []
        with open("xe.mcf") as f:
            for line in f:
                mcf_data.append(line.strip().split())

        for idx, line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == "Atom_Info":
                    atom_section_start = idx

        # Check a some atom info
        assert mcf_data[atom_section_start + 1][0] == "1"
        assert mcf_data[atom_section_start + 2][1] == "Xe"
        assert mcf_data[atom_section_start + 2][2] == "Xe"
        assert mcf_data[atom_section_start + 2][5] == "Mie"
        assert np.isclose(
            float(mcf_data[atom_section_start + 2][3]),
            top.sites[0].mass.in_units(u.amu).value,
        )
        assert np.isclose(
            float(mcf_data[atom_section_start + 2][4]),
            top.sites[0].charge.in_units(u.elementary_charge).value,
        )
        assert np.isclose(
            float(mcf_data[atom_section_start + 2][6]),
            (top.sites[0].atom_type.parameters["epsilon"] / u.kb)
            .in_units(u.K)
            .value,
        )
        assert np.isclose(
            float(mcf_data[atom_section_start + 2][7]),
            top.sites[0]
            .atom_type.parameters["sigma"]
            .in_units(u.Angstrom)
            .value,
        )
        assert np.isclose(
            float(mcf_data[atom_section_start + 2][8]),
            top.sites[0].atom_type.parameters["n"],
        )
        assert np.isclose(
            float(mcf_data[atom_section_start + 2][9]),
            top.sites[0].atom_type.parameters["m"],
        )

    def test_modified_potentials(self, n_typed_ar_system):
        top = n_typed_ar_system(n_sites=1)

        top.atom_types[0].set_expression("sigma + epsilon*r")

        with pytest.raises(EngineIncompatibilityError):
            top.save("out.mcf")

        alternate_lj = "4*epsilon*sigma**12/r**12 - 4*epsilon*sigma**6/r**6"
        top.atom_types[0].set_expression(alternate_lj)

        top.save("ar.mcf")

    def test_scaling_factors(self, n_typed_ar_system):
        top = n_typed_ar_system()
        top.save("ar.mcf")
        mcf_data = []
        with open("ar.mcf") as f:
            for line in f:
                mcf_data.append(line.strip().split())
        assert np.allclose(float(mcf_data[-5][0]), 0.0)
        assert np.allclose(float(mcf_data[-5][1]), 0.0)
        assert np.allclose(float(mcf_data[-5][2]), 0.5)
        assert np.allclose(float(mcf_data[-5][3]), 1.0)
        assert np.allclose(float(mcf_data[-4][0]), 0.0)
        assert np.allclose(float(mcf_data[-4][1]), 0.0)
        assert np.allclose(float(mcf_data[-4][2]), 0.5)
        assert np.allclose(float(mcf_data[-4][3]), 1.0)

        top.scaling_factors = {
            "nonBonded12Scale": 0.1,
            "nonBonded13Scale": 0.2,
            "nonBonded14Scale": 0.5,
            "electrostatics12Scale": 0.2,
            "electrostatics13Scale": 0.4,
            "electrostatics14Scale": 0.6,
        }

        top.save("ar.mcf", overwrite=True)
        mcf_data = []
        with open("ar.mcf") as f:
            for line in f:
                mcf_data.append(line.strip().split())
        assert np.allclose(float(mcf_data[-5][0]), 0.1)
        assert np.allclose(float(mcf_data[-5][1]), 0.2)
        assert np.allclose(float(mcf_data[-5][2]), 0.5)
        assert np.allclose(float(mcf_data[-5][3]), 1.0)
        assert np.allclose(float(mcf_data[-4][0]), 0.2)
        assert np.allclose(float(mcf_data[-4][1]), 0.4)
        assert np.allclose(float(mcf_data[-4][2]), 0.6)
        assert np.allclose(float(mcf_data[-4][3]), 1.0)
