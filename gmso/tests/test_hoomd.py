import forcefield_utilities as ffutils
import hoomd
import numpy as np
import pytest
import unyt as u
from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield

from gmso import ForceField
from gmso.external import from_mbuild
from gmso.external.convert_hoomd import to_hoomd_forcefield, to_hoomd_snapshot
from gmso.parameterization import apply
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path
from gmso.utils.io import has_hoomd, has_mbuild, import_
from gmso.utils.sorting import sort_connection_strings

if has_hoomd:
    hoomd = import_("hoomd")
    hoomd_version = hoomd.version.version.split(".")

if has_mbuild:
    mb = import_("mbuild")


def run_hoomd_nvt(snapshot, forces, vhoomd=4):
    cpu = hoomd.device.CPU()
    sim = hoomd.Simulation(device=cpu)
    sim.create_state_from_snapshot(snapshot)

    integrator = hoomd.md.Integrator(dt=0.001)
    integrator.forces = list(set().union(*forces.values()))

    temp = 300 * u.K
    kT = temp.to_equivalent("kJ/mol", "thermal").value
    if vhoomd == 4:
        thermostat = hoomd.md.methods.thermostats.MTTK(kT=kT, tau=1.0)
        nvt = hoomd.md.methods.ConstantVolume(
            thermostat=thermostat, filter=hoomd.filter.All()
        )
    elif vhoomd == 3:
        nvt = hoomd.md.methods.NVT(kT=kT, filter=hoomd.filter.All(), tau=1.0)
    else:
        raise ImportError("Wrong version of hoomd.")
    integrator.methods.append(nvt)
    sim.operations.integrator = integrator

    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT)
    thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(
        filter=hoomd.filter.All()
    )

    sim.operations.computes.append(thermodynamic_properties)
    return sim


@pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
@pytest.mark.skipif(not has_mbuild, reason="mbuild not installed")
class TestGsd(BaseTest):
    def test_mbuild_comparison(self, oplsaa_forcefield):
        compound = mb.load("CCC", smiles=True)
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=2)
        base_units = {
            "mass": u.g / u.mol,
            "length": u.nm,
            "energy": u.kJ / u.mol,
        }

        top = from_mbuild(com_box)
        top.identify_connections()
        top = apply(top, oplsaa_forcefield, remove_untyped=True)

        gmso_snapshot, _ = to_hoomd_snapshot(top, base_units=base_units)
        gmso_forces, _ = to_hoomd_forcefield(
            top,
            r_cut=1.4,
            base_units=base_units,
            pppm_kwargs={"resolution": (64, 64, 64), "order": 7},
        )

        integrator_forces = list()
        for cat in gmso_forces:
            for force in gmso_forces[cat]:
                integrator_forces.append(force)

        import foyer

        oplsaa = foyer.Forcefield(name="oplsaa")
        structure = oplsaa.apply(com_box)

        d = 10
        e = 1 / 4.184
        m = 0.9999938574

        mb_snapshot, mb_forcefield, _ = create_hoomd_forcefield(
            structure,
            ref_distance=d,
            ref_energy=e,
            ref_mass=m,
            r_cut=1.4,
            init_snap=None,
            pppm_kwargs={"Nx": 64, "Ny": 64, "Nz": 64, "order": 7},
        )

        assert mb_snapshot.particles.N == gmso_snapshot.particles.N
        assert np.allclose(
            mb_snapshot.particles.position, gmso_snapshot.particles.position
        )
        assert mb_snapshot.bonds.N == gmso_snapshot.bonds.N
        assert mb_snapshot.angles.N == gmso_snapshot.angles.N
        assert mb_snapshot.dihedrals.N == gmso_snapshot.dihedrals.N

        sorted_gmso_ff = sorted(
            integrator_forces, key=lambda cls: str(cls.__class__)
        )
        sorted_mbuild_ff = sorted(
            mb_forcefield, key=lambda cls: str(cls.__class__)
        )

        for mb_force, gmso_force in zip(sorted_mbuild_ff, sorted_gmso_ff):
            if isinstance(
                mb_force,
                (
                    hoomd.md.long_range.pppm.Coulomb,
                    hoomd.md.pair.pair.LJ,
                    hoomd.md.special_pair.LJ,
                    hoomd.md.pair.pair.Ewald,
                    hoomd.md.special_pair.Coulomb,
                ),
            ):
                continue
            keys = mb_force.params.param_dict.keys()
            for key in keys:
                gmso_key = key.replace("opls_135", "CT")
                gmso_key = gmso_key.replace("opls_136", "CT")
                gmso_key = gmso_key.replace("opls_140", "HC")
                gmso_key = "-".join(
                    sort_connection_strings(gmso_key.split("-"))
                )
                mb_params = mb_force.params.param_dict[key]
                gmso_params = gmso_force.params.param_dict[gmso_key]
                variables = mb_params.keys()
                for var in variables:
                    assert np.isclose(mb_params[var], gmso_params[var])

    @pytest.mark.skipif(
        int(hoomd_version[0]) < 4, reason="Unsupported features in HOOMD 3"
    )
    def test_hoomd4_simulation(self):
        compound = mb.load("CCC", smiles=True)
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=2)
        base_units = {
            "mass": u.g / u.mol,
            "length": u.nm,
            "energy": u.kJ / u.mol,
        }

        top = from_mbuild(com_box)
        top.identify_connections()
        oplsaa = ffutils.FoyerFFs().load("oplsaa").to_gmso_ff()
        top = apply(top, oplsaa, remove_untyped=True)

        gmso_snapshot, _ = to_hoomd_snapshot(top, base_units=base_units)
        gmso_forces, _ = to_hoomd_forcefield(
            top,
            r_cut=1.4,
            base_units=base_units,
            pppm_kwargs={"resolution": (64, 64, 64), "order": 7},
        )

        sim = run_hoomd_nvt(gmso_snapshot, gmso_forces)
        sim.run(100)

    @pytest.mark.skipif(
        int(hoomd_version[0]) < 4, reason="Deprecated features in HOOMD 4"
    )
    def test_hoomd4_simulation_auto_scaled(self):
        compound = mb.load("CCC", smiles=True)
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=2)
        base_units = {
            "mass": u.g / u.mol,
            "length": u.nm,
            "energy": u.kJ / u.mol,
        }

        top = from_mbuild(com_box)
        top.identify_connections()
        oplsaa = ffutils.FoyerFFs().load("oplsaa").to_gmso_ff()
        top = apply(top, oplsaa, remove_untyped=True)

        gmso_snapshot, _ = to_hoomd_snapshot(
            top,
            base_units=base_units,
            auto_scale=True,
        )
        gmso_forces, _ = to_hoomd_forcefield(
            top,
            r_cut=1.4,
            base_units=base_units,
            pppm_kwargs={"resolution": (64, 64, 64), "order": 7},
            auto_scale=True,
        )

        sim = run_hoomd_nvt(gmso_snapshot, gmso_forces)
        sim.run(100)

    @pytest.mark.skipif(
        int(hoomd_version[0]) >= 4, reason="Deprecated features in HOOMD 4"
    )
    def test_hoomd3_simulation(self):
        compound = mb.load("CCC", smiles=True)
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=2)
        base_units = {
            "mass": u.g / u.mol,
            "length": u.nm,
            "energy": u.kJ / u.mol,
        }

        top = from_mbuild(com_box)
        top.identify_connections()
        oplsaa = ffutils.FoyerFFs().load("oplsaa").to_gmso_ff()
        top = apply(top, oplsaa, remove_untyped=True)

        gmso_snapshot, snapshot_base_units = to_hoomd_snapshot(
            top, base_units=base_units
        )
        gmso_forces, forces_base_units = to_hoomd_forcefield(
            top,
            r_cut=1.4,
            base_units=base_units,
            pppm_kwargs={"resolution": (64, 64, 64), "order": 7},
        )

        sim = run_hoomd_nvt(gmso_snapshot, gmso_forces, vhoomd=3)
        sim.run(100)

    @pytest.mark.skipif(
        int(hoomd_version[0]) >= 4, reason="Deprecated features in HOOMD 4"
    )
    def test_hoomd3_simulation_auto_scaled(self):
        compound = mb.load("CCC", smiles=True)
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=2)
        base_units = {
            "mass": u.g / u.mol,
            "length": u.nm,
            "energy": u.kJ / u.mol,
        }

        top = from_mbuild(com_box)
        top.identify_connections()
        oplsaa = ffutils.FoyerFFs().load("oplsaa").to_gmso_ff()
        top = apply(top, oplsaa, remove_untyped=True)

        gmso_snapshot, snapshot_base_units = to_hoomd_snapshot(
            top,
            base_units=base_units,
            auto_scale=True,
        )
        gmso_forces, forces_base_units = to_hoomd_forcefield(
            top,
            r_cut=1.4,
            base_units=base_units,
            pppm_kwargs={"resolution": (64, 64, 64), "order": 7},
            auto_scale=True,
        )

        integrator_forces = list()
        for cat in gmso_forces:
            for force in gmso_forces[cat]:
                integrator_forces.append(force)

        temp = 300 * u.K
        kT = temp.to_equivalent("kJ/mol", "thermal").value

        cpu = hoomd.device.CPU()
        sim = hoomd.Simulation(device=cpu)
        sim.create_state_from_snapshot(gmso_snapshot)

        integrator = hoomd.md.Integrator(dt=0.001)
        # cell = hoomd.md.nlist.Cell(buffer=0.4)
        integrator.forces = integrator_forces
        # integrator.forces = mb_forcefield

        nvt = hoomd.md.methods.NVT(kT=kT, filter=hoomd.filter.All(), tau=1.0)
        integrator.methods.append(nvt)
        sim.operations.integrator = integrator

        sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT)
        thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(
            filter=hoomd.filter.All()
        )

        sim.operations.computes.append(thermodynamic_properties)
        sim.run(100)

    def test_diff_base_units(self):
        compound = mb.load("CC", smiles=True)
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=2)
        base_units = {
            "mass": u.amu,
            "length": u.nm,
            "energy": u.kJ / u.mol,
        }

        top = from_mbuild(com_box)
        top.identify_connections()
        oplsaa = ffutils.FoyerFFs().load("oplsaa").to_gmso_ff()
        top = apply(top, oplsaa, remove_untyped=True)

        gmso_snapshot, snapshot_base_units = to_hoomd_snapshot(
            top, base_units=base_units
        )
        gmso_forces, forces_base_units = to_hoomd_forcefield(
            top,
            r_cut=1.4,
            base_units=base_units,
            pppm_kwargs={"resolution": (64, 64, 64), "order": 7},
        )

    def test_default_units(self):
        compound = mb.load("CC", smiles=True)
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=2)
        base_units = {
            "mass": u.amu,
            "length": u.nm,
            "energy": u.kJ / u.mol,
        }

        top = from_mbuild(com_box)
        top.identify_connections()
        oplsaa = ffutils.FoyerFFs().load("oplsaa").to_gmso_ff()
        top = apply(top, oplsaa, remove_untyped=True)

        gmso_snapshot, snapshot_base_units = to_hoomd_snapshot(top)
        gmso_forces, forces_base_units = to_hoomd_forcefield(
            top=top,
            r_cut=1.4,
            pppm_kwargs={"resolution": (64, 64, 64), "order": 7},
        )

    def test_ff_zero_parameter(self):
        ethane = mb.lib.molecules.Ethane()
        top = from_mbuild(ethane)
        top.identify_connections()
        ff_zero_param = (
            ffutils.FoyerFFs()
            .load(get_path("ethane_zero_parameter.xml"))
            .to_gmso_ff()
        )
        top = apply(top, ff_zero_param, remove_untyped=True)
        base_units = {
            "mass": u.g / u.mol,
            "length": u.nm,
            "energy": u.kJ / u.mol,
        }
        gmso_forces, forces_base_units = to_hoomd_forcefield(
            top,
            r_cut=1.4,
            base_units=base_units,
            pppm_kwargs={"resolution": (64, 64, 64), "order": 7},
        )
        integrator_forces = list()
        for cat in gmso_forces:
            for force in gmso_forces[cat]:
                integrator_forces.append(force)
        for force in integrator_forces:
            if isinstance(force, hoomd.md.pair.LJ):
                keys = force.params.param_dict.keys()
                for key in keys:
                    if "opls_135" in list(key):
                        params = force.params.param_dict[key]
                        variables = params.keys()
                        for var in variables:
                            assert params[var] == 0.0

    def test_zero_charges(self):
        compound = mb.load("CC", smiles=True)
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=2)
        base_units = {
            "mass": u.amu,
            "length": u.nm,
            "energy": u.kJ / u.mol,
        }

        top = from_mbuild(com_box)
        top.identify_connections()
        oplsaa = ffutils.FoyerFFs().load("oplsaa").to_gmso_ff()
        top = apply(top, oplsaa, remove_untyped=True)
        for site in top.sites:
            site.charge = 0

        gmso_forces, _ = to_hoomd_forcefield(
            top=top,
            r_cut=1.4,
        )
        for cat in gmso_forces:
            for force in gmso_forces[cat]:
                assert not isinstance(force, hoomd.md.pair.pair.Ewald)
                assert not isinstance(force, hoomd.md.long_range.pppm.Coulomb)
                assert not isinstance(force, hoomd.md.special_pair.Coulomb)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    @pytest.mark.skipif(not has_mbuild, reason="mbuild not installed")
    @pytest.mark.skipif(
        int(hoomd_version[0]) <= 3.8, reason="Deprecated features in HOOMD 4"
    )
    def test_gaff_sim(self, gaff_forcefield):
        base_units = {
            "mass": u.g / u.mol,
            "length": u.nm,
            "energy": u.kJ / u.mol,
        }
        benzene = mb.load("c1ccccc1", smiles=True)
        benzene.box = mb.Box([5, 5, 5])
        top = benzene.to_gmso()
        parameterized_top = apply(
            top, gaff_forcefield, identify_connections=True
        )
        assert parameterized_top.is_fully_typed

        snap, _ = to_hoomd_snapshot(parameterized_top, base_units)
        forces, _ = to_hoomd_forcefield(
            parameterized_top,
            r_cut=1.4,
            base_units=base_units,
            pppm_kwargs={"resolution": (64, 64, 64), "order": 7},
        )
        assert forces
        assert snap

        sim = run_hoomd_nvt(snap, forces)
        sim.run(100)

    def test_forces_connections_match(self):
        compound = mb.load("CC", smiles=True)
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=1)
        base_units = {
            "mass": u.amu,
            "length": u.nm,
            "energy": u.kJ / u.mol,
        }
        top = com_box.to_gmso()
        top.identify_connections()
        ethaneFF = ForceField(get_path("alkanes.xml"))

        top = apply(top, ethaneFF, remove_untyped=True)

        snapshot, _ = to_hoomd_snapshot(top, base_units=base_units)
        assert "CT-HC" in snapshot.bonds.types

        forces, _ = to_hoomd_forcefield(
            top=top, r_cut=1.4, base_units=base_units
        )
        assert "CT-HC" in forces["bonds"][0].params.keys()

    def test_forces_wildcards(self):
        compound = mb.load("CCCC", smiles=True)
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=1)
        base_units = {
            "mass": u.amu,
            "length": u.nm,
            "energy": u.kJ / u.mol,
        }
        top = com_box.to_gmso()
        top.identify_connections()
        ethaneFF = ForceField(get_path("alkanes_wildcards.xml"))
        top = apply(top, ethaneFF, remove_untyped=True)

        snapshot, _ = to_hoomd_snapshot(top, base_units=base_units)
        assert "CT-HC" in snapshot.bonds.types

        forces, _ = to_hoomd_forcefield(
            top=top, r_cut=1.4, base_units=base_units
        )
        assert "CT-CT-CT-HC" in list(forces["dihedrals"][0].params)
        for conntype in snapshot.dihedrals.types:
            assert conntype in list(forces["dihedrals"][0].params)
