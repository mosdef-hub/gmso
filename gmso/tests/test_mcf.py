import re
import subprocess

import mbuild as mb
import numpy as np
import pytest
import unyt as u

from gmso.core.forcefield import ForceField
from gmso.exceptions import EngineIncompatibilityError
from gmso.external import from_mbuild
from gmso.formats.mcf import write_mcf
from gmso.parameterization import apply
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path
from gmso.utils.conversions import convert_ryckaert_to_fourier
from gmso.utils.io import get_fn, has_cassandra, has_parmed, import_

if has_cassandra:
    mc = import_("mosdef_cassandra")

if has_parmed:
    pmd = import_("parmed")


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
        if len(line) > 1:
            if line[1] == "Intra_Scaling":
                mcf_idx["Intra_Scaling"] = idx

    return mcf_data, mcf_idx


def run_cassandra(cassandra, inp_file):
    """Calls Cassandra. Taken from mosdef_cassandra v0.3.2"""
    cassandra_cmd = "{cassandra} {inp_file}".format(
        cassandra=cassandra, inp_file=inp_file
    )
    p = subprocess.Popen(
        cassandra_cmd,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    out, err = p.communicate()

    if p.returncode != 0 or "error" in err.lower() or "error" in out.lower():
        return 1, out, err
    return 0, out, err


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
            top.sites[0].atom_type.mass.in_units(u.amu).value,
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
            top.sites[0].atom_type.mass.in_units(u.amu).value,
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

    def test_write_single_fragment_two_atoms(self):
        """
        The main purpose of his test is to check for valid
        fragment information ouput for molecules that have
        no angles.
        """
        ethane = mb.load(get_fn("ethane_ua.mol2"))
        top = from_mbuild(ethane)
        ff = ForceField(get_path("ethane-rigid.xml"))
        top.identify_connections()
        apply(top, ff, remove_untyped=True)
        write_mcf(top, "ethane-rigid.mcf")

        mcf_data, mcf_idx = parse_mcf("ethane-rigid.mcf")

        assert is_charge_neutral(mcf_data, mcf_idx)

        # Assert number of fragments
        assert mcf_data[mcf_idx["Fragment_Info"] + 1][0] == "1"
        # Assert number of atoms in the first fragment
        assert mcf_data[mcf_idx["Fragment_Info"] + 2][1] == "2"
        # Assert atom IDs in the first fragment
        assert mcf_data[mcf_idx["Fragment_Info"] + 2][2] == "1"
        assert mcf_data[mcf_idx["Fragment_Info"] + 2][3] == "2"
        assert mcf_data[mcf_idx["Fragment_Connectivity"] + 1][0] == "0"

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

    def test_typed_ethylene(self):
        ethylene = mb.load("C=C", smiles=True)
        top = from_mbuild(ethylene)
        ff = ForceField(get_path("ethylene.xml"))
        top.identify_connections()
        apply(top, ff, remove_untyped=True)
        write_mcf(top, "ethylene.mcf")

        mcf_data, mcf_idx = parse_mcf("ethylene.mcf")

        assert is_charge_neutral(mcf_data, mcf_idx)

        # Check atom info
        assert mcf_data[mcf_idx["Atom_Info"] + 1][0] == "6"
        assert mcf_data[mcf_idx["Atom_Info"] + 2][1] == "opls_143"
        assert mcf_data[mcf_idx["Atom_Info"] + 2][2] == "C"
        assert mcf_data[mcf_idx["Atom_Info"] + 2][5] == "LJ"
        assert np.isclose(
            float(mcf_data[mcf_idx["Atom_Info"] + 2][3]),
            top.sites[0].atom_type.mass.in_units(u.amu).value,
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
        assert mcf_data[mcf_idx["Bond_Info"] + 1][0] == "5"
        assert mcf_data[mcf_idx["Bond_Info"] + 2][3] == "fixed"
        assert np.isclose(
            float(mcf_data[mcf_idx["Bond_Info"] + 2][4]),
            top.bonds[0]
            .bond_type.parameters["r_eq"]
            .in_units(u.Angstrom)
            .value,
        )

        # Check angle section
        assert mcf_data[mcf_idx["Angle_Info"] + 1][0] == "6"
        assert mcf_data[mcf_idx["Angle_Info"] + 2][4] == "harmonic"
        assert np.isclose(
            float(mcf_data[mcf_idx["Angle_Info"] + 2][6]),
            top.angles[0]
            .angle_type.parameters["theta_eq"]
            .in_units(u.degree)
            .value,
        )

        assert np.isclose(
            2.0 * float(mcf_data[mcf_idx["Angle_Info"] + 2][5]),
            (top.angles[0].angle_type.parameters["k"] / u.kb)
            .in_units(u.K / u.radian**2)
            .value,
        )

        # Check dihedral section
        assert mcf_data[mcf_idx["Dihedral_Info"] + 1][0] == "4"
        dihedral_style = mcf_data[mcf_idx["Dihedral_Info"] + 2][5].lower()
        assert dihedral_style.lower() == "opls"

        k = list(ff.dihedral_types.keys())
        dihedral_type = ff.dihedral_types[k[0]]
        ff_coeffs = convert_ryckaert_to_fourier(
            dihedral_type
        ).parameters.values()

        # We need to drop the last coefficient in the fourier GMSO dihedral type
        # (the term that has cos(4*phi) since it is not supported in the equivalent OPLS
        # dihedral type in Cassandra)

        ff_coeffs = np.array([float(x) for x in ff_coeffs])[:-1]
        mcf_coeffs = np.array(
            [float(x) for x in mcf_data[mcf_idx["Dihedral_Info"] + 2][6:]]
        )

        assert np.allclose(ff_coeffs, 2.0 * mcf_coeffs)

    def test_fixed_angles(self, typed_tip3p_rigid_system):
        top = typed_tip3p_rigid_system
        write_mcf(top, "tip3p-rigid.mcf")

        mcf_data, mcf_idx = parse_mcf("tip3p-rigid.mcf")

        assert is_charge_neutral(mcf_data, mcf_idx)

        # Check angle section
        assert mcf_data[mcf_idx["Angle_Info"] + 1][0] == "1"
        assert mcf_data[mcf_idx["Angle_Info"] + 2][4] == "fixed"
        assert np.isclose(
            float(mcf_data[mcf_idx["Angle_Info"] + 2][5]),
            top.angles[0]
            .angle_type.parameters["theta_eq"]
            .in_units(u.degree)
            .value,
        )

    def test_top_with_n_ar_system(self, n_typed_ar_system):
        top = n_typed_ar_system(n_sites=10)
        write_mcf(top, "top-10ar.mcf")

        mcf_data, mcf_idx = parse_mcf("top-10ar.mcf")

        assert is_charge_neutral(mcf_data, mcf_idx)

        assert mcf_data[mcf_idx["Atom_Info"] + 1][0] == "1"
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

    def test_top_with_mixture(self):
        pass

    @pytest.mark.skipif(not has_cassandra, reason="cassandra is not installed")
    def test_in_cassandra(self, typed_ethane):
        """
        This test runs a single point energy calculation in Cassandra using an MCF
        generated by gmso and compare the total energy
        to the energy of a simulation run with a MCF file generated using mosdef_cassandra
        (which involves using a parmed.Structure)
        """
        from mosdef_cassandra.utils.detect import detect_cassandra_binaries
        from mosdef_cassandra.writers.writers import write_input

        from gmso.external.convert_parmed import to_parmed

        # First run the mosdef_cassandra simulation. Mosdef_cassandra generates an input file
        # as well as an MCF. Later, we will use a mosdef_cassandra
        # input file as a template and replace the MCF file with the GMSO MCF file

        box = mb.Box([3.0, 3.0, 3.0])
        species = to_parmed(typed_ethane)
        system = mc.System([box], [species], mols_to_add=[[5]])
        ensemble = "nvt"
        moveset = mc.MoveSet(ensemble, [species])
        mc.run(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=0,
            temperature=300.0 * u.K,
            run_name="nvt_mbuild",
            seeds=[12356, 64321],
        )

        py, fraglib_setup, cassandra = detect_cassandra_binaries()

        # TODO: not sure why the cassandra MCF writer of mBuild
        # outputs a different intramolecular exclusions relative
        # to the GMSO writer. This is a temporary fix to make the
        # test pass. We should investigate this further.
        # Also, try to use the function top.set_lj_scale() and see how
        # to update subtopologies.

        write_mcf(typed_ethane, "gmso.mcf")
        mcf_data, mcf_idx = parse_mcf("gmso.mcf")
        mcf_data[mcf_idx["Intra_Scaling"] + 1][0:4] = [
            "0.0",
            "0.0",
            "0.0",
            "1.0",
        ]
        mcf_data[mcf_idx["Intra_Scaling"] + 2][0:4] = [
            "0.0",
            "0.0",
            "0.0",
            "1.0",
        ]
        with open("gmso.mcf", mode="w") as f:
            for line in mcf_data:
                f.write(" ".join(line) + "\n")

        inp_file = write_input(
            system=system,
            moveset=moveset,
            run_type="equilibration",
            run_length=0,
            temperature=300.0 * u.K,
            run_name="nvt_gmso",
            seeds=[12356, 64321],
        )

        with open(inp_file, mode="r") as f:
            lines = f.read()
            lines = lines.replace("species1.mcf", "gmso.mcf")
            # The fragment files section is empty unless
            # restart option is used in mosdef_cassandra.
            # See the mosdef_cassandra.writers.inp_functions.generate_input
            lines = lines.replace(
                "# Fragment_Files",
                "# Fragment_Files\nspecies1/frag1/frag1.dat  1\nspecies1/frag2/frag2.dat  2\n",
            )

        with open(inp_file, mode="w") as f:
            f.writelines(lines)

        # Run the simulation with the GMSO MCF file
        code, out, err = run_cassandra(cassandra, inp_file)

        assert code == 0
        assert "complete" in out

        # Parse log files
        with open("nvt_gmso.out.log", mode="r") as f:
            lines = f.readlines()
            energy = None
            for line in lines:
                if "Total system energy" in line:
                    energy = float(line.split()[-1])
                    break

        with open("nvt_mbuild.out.log", mode="r") as f:
            lines = f.readlines()
            energy_ref = 0.0
            for line in lines:
                if "Total system energy" in line:
                    energy_ref = float(line.split()[-1])
                    break

        assert np.isclose(energy, energy_ref, rtol=1e-3)

    @pytest.mark.skipif(not has_parmed, reason="ParmEd is not installed")
    def test_parmed_vs_gmso(self, parmed_ethane):
        """
        This test compares the output of a MCF file generated
        by gmso to the output of a MCF file generated by parmed.
        The Parmed MCF file is generated through mbuild.
        """
        from gmso.external.convert_parmed import from_parmed

        mb.formats.cassandramcf.write_mcf(
            parmed_ethane,
            "parmed-ethane.mcf",
            dihedral_style="opls",
            angle_style="harmonic",
        )

        top = from_parmed(parmed_ethane)
        write_mcf(top, "gmso-ethane.mcf")

        mcf_data_pmd, mcf_idx_pmd = parse_mcf("parmed-ethane.mcf")
        mcf_data_gmso, mcf_idx_gmso = parse_mcf("gmso-ethane.mcf")
        skip_lines = [3]
        float_pattern = r"[+-]?[0-9]*[.][0-9]*"
        for i, (line_pmd, line_gmso) in enumerate(
            zip(mcf_data_pmd, mcf_data_gmso)
        ):
            if i in skip_lines:
                continue

            pmd_parms = np.array(
                re.findall(float_pattern, " ".join(line_pmd)), dtype=np.float64
            )
            gmso_parms = np.array(
                re.findall(float_pattern, " ".join(line_gmso)), dtype=np.float64
            )

            assert np.allclose(pmd_parms, gmso_parms, rtol=1e-3)

        # TODO: can we read top.mcf into the simulation somehow?
        # TODO: or use it to validate that the simulation was done
        # correctly?

    def test_untyped_top(self):
        pass

    def test_top_with_ring(self, typed_benzene_ua_system):
        top = typed_benzene_ua_system
        write_mcf(top, "benzene-ua.mcf")

        mcf_data, mcf_idx = parse_mcf("benzene-ua.mcf")

        assert is_charge_neutral(mcf_data, mcf_idx)

        assert mcf_data[mcf_idx["Atom_Info"] + 1][0] == "6"

        for idx in range(0, 6):
            last_label = mcf_data[mcf_idx["Atom_Info"] + 2 + idx][-1]
            assert last_label == "ring"

        assert mcf_data[mcf_idx["Bond_Info"] + 1][0] == "6"
        assert mcf_data[mcf_idx["Angle_Info"] + 1][0] == "6"
        assert mcf_data[mcf_idx["Dihedral_Info"] + 1][0] == "6"

        assert mcf_data[mcf_idx["Fragment_Info"] + 1][0] == "1"
        frag_atoms = mcf_data[mcf_idx["Fragment_Info"] + 2][1:]
        assert set(frag_atoms) == set([str(i) for i in range(1, 7)])
