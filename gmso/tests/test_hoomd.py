import forcefield_utilities as ffutils
import hoomd
import numpy as np
import pytest
import unyt as u

from gmso import ForceField
from gmso.external import from_mbuild
from gmso.external.convert_hoomd import (
    to_gsd_snapshot,
    to_hoomd_forcefield,
    to_hoomd_snapshot,
)
from gmso.parameterization import apply
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path
from gmso.utils.io import has_hoomd, has_mbuild, import_

if has_hoomd:
    hoomd = import_("hoomd")
    hoomd_version = hoomd.version.version.split(".")

if has_mbuild:
    mb = import_("mbuild")


def run_hoomd_nvt(snapshot, forces, vhoomd=4):
    cpu = hoomd.device.CPU()
    sim = hoomd.Simulation(device=cpu)
    sim.create_state_from_snapshot(snapshot)

    integrator = hoomd.md.Integrator(dt=0.0001)
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


@pytest.mark.skipif(not has_mbuild, reason="mbuild not installed")
class TestGsd(BaseTest):
    def test_rigid_bodies(self):
        ethane = mb.lib.molecules.Ethane()
        box = mb.fill_box(ethane, n_compounds=2, box=[2, 2, 2])
        top = from_mbuild(box)
        for site in top.sites:
            site.molecule.isrigid = True
        top.identify_connections()

        top_no_rigid = from_mbuild(box)
        top_no_rigid.identify_connections()

        rigid_ids = [site.molecule.number for site in top.sites]
        assert set(rigid_ids) == {0, 1}

        snapshot, refs, rigid = to_gsd_snapshot(top)
        snapshot.validate()
        snapshot_no_rigid, refs, _ = to_gsd_snapshot(top_no_rigid)
        # Check that snapshot has rigid particles added
        assert "Ethane" in snapshot.particles.types
        assert "Ethane" not in snapshot_no_rigid.particles.types
        assert snapshot.particles.N - 2 == snapshot_no_rigid.particles.N
        assert np.array_equal(snapshot.particles.typeid[:2], np.array([0, 0]))
        assert np.array_equal(
            snapshot.particles.body[2:10], np.array([0] * ethane.n_particles)
        )
        assert np.array_equal(
            snapshot.particles.body[10:], np.array([1] * ethane.n_particles)
        )
        assert np.allclose(snapshot.particles.mass[0], ethane.mass, atol=1e-2)
        # Check that topology isn't changed by using rigid bodies
        assert snapshot.bonds.N == snapshot_no_rigid.bonds.N
        assert snapshot.angles.N == snapshot_no_rigid.angles.N
        assert snapshot.dihedrals.N == snapshot_no_rigid.dihedrals.N
        # Check if particle indices in snapshot groups are adjusted by N rigid bodies
        for group1, group2 in zip(snapshot.bonds.group, snapshot_no_rigid.bonds.group):
            assert np.array_equal(np.array(group1), np.array(group2) + 2)
        for group1, group2 in zip(
            snapshot.angles.group, snapshot_no_rigid.angles.group
        ):
            assert np.array_equal(np.array(group1), np.array(group2) + 2)
        for group1, group2 in zip(
            snapshot.dihedrals.group, snapshot_no_rigid.dihedrals.group
        ):
            assert np.array_equal(np.array(group1), np.array(group2) + 2)
        for group1, group2 in zip(
            snapshot.impropers.group, snapshot_no_rigid.impropers.group
        ):
            assert np.array_equal(np.array(group1), np.array(group2) + 2)

    def test_multiple_rigid_bodies(self, gaff_forcefield):
        ethane = mb.lib.molecules.Ethane()
        benzene = mb.load("c1ccccc1", smiles=True)
        benzene.name = "Benzene"
        box = mb.fill_box([ethane, benzene], n_compounds=[1, 1], box=[2, 2, 2])
        top = from_mbuild(box)
        for site in top.sites:
            site.molecule.isrigid = True
        apply(top, gaff_forcefield, identify_connections=True)

        rigid_ids = [site.molecule.number for site in top.sites]
        assert set(rigid_ids) == {0, 1}

        snapshot, refs, rigid = to_gsd_snapshot(top)
        snapshot.validate()
        ethane_types = list(
            set(
                [
                    site.atom_type.name
                    for site in top.sites
                    if site.molecule.name == "Ethane"
                ]
            )
        )
        benzene_types = list(
            set(
                [
                    site.atom_type.name
                    for site in top.sites
                    if site.molecule.name == "Benzene"
                ]
            )
        )
        for _type in ethane_types:
            assert _type in rigid.body["Ethane"]["constituent_types"]

        for _type in benzene_types:
            assert _type in rigid.body["Benzene"]["constituent_types"]

        assert np.array_equal(np.zeros(8), snapshot.particles.body[2:10])
        assert np.array_equal(np.ones(12), snapshot.particles.body[10:])
        assert round(float(snapshot.particles.mass[0]), 4) == round(
            float(ethane.mass), 4
        )
        assert round(float(snapshot.particles.mass[1]), 4) == round(
            float(benzene.mass), 4
        )

    def test_rigid_bodies_mix(self):
        ethane = mb.lib.molecules.Ethane()
        ethane.name = "ethane"
        benzene = mb.load("c1ccccc1", smiles=True)
        benzene.name = "benzene"
        box = mb.fill_box([benzene, ethane], n_compounds=[1, 1], box=[2, 2, 2])
        top = from_mbuild(box)
        for site in top.sites:
            if site.molecule.name == "benzene":
                site.molecule.isrigid = True

        snapshot, refs, rigid = to_gsd_snapshot(top)
        snapshot.validate()
        assert snapshot.particles.typeid[0] == 0
        assert snapshot.particles.N == top.n_sites + 1
        assert np.array_equal(
            snapshot.particles.body[-ethane.n_particles :],
            np.array([-1] * ethane.n_particles),
        )
        assert np.array_equal(
            snapshot.particles.body[1 : benzene.n_particles + 1],
            np.array([0] * benzene.n_particles),
        )

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

        gmso_snapshot, _, _ = to_hoomd_snapshot(top, base_units=base_units)
        gmso_snapshot.validate()
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

        gmso_snapshot, _, _ = to_hoomd_snapshot(
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

        gmso_snapshot, snapshot_base_units, _ = to_hoomd_snapshot(
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

        gmso_snapshot, snapshot_base_units, _ = to_hoomd_snapshot(
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

        gmso_snapshot, snapshot_base_units, _ = to_hoomd_snapshot(
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

        gmso_snapshot, snapshot_base_units, _ = to_hoomd_snapshot(top)
        gmso_forces, forces_base_units = to_hoomd_forcefield(
            top=top,
            r_cut=1.4,
            pppm_kwargs={"resolution": (64, 64, 64), "order": 7},
            base_units=base_units,
        )

    def test_ff_zero_parameter(self):
        ethane = mb.lib.molecules.Ethane()
        top = from_mbuild(ethane)
        top.identify_connections()
        ff_zero_param = (
            ffutils.FoyerFFs().load(get_path("ethane_zero_parameter.xml")).to_gmso_ff()
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
        int(hoomd_version[0]) < 4.5, reason="No periodic impropers in hoomd < 4.5"
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
        parameterized_top = apply(top, gaff_forcefield, identify_connections=True)
        assert parameterized_top.is_fully_typed

        snap, _, _ = to_hoomd_snapshot(parameterized_top, base_units)
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

        snapshot, _, _ = to_hoomd_snapshot(top, base_units=base_units)
        assert "CT-HC" in snapshot.bonds.types

        forces, _ = to_hoomd_forcefield(top=top, r_cut=1.4, base_units=base_units)
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

        snapshot, _, _ = to_hoomd_snapshot(top, base_units=base_units)
        assert "CT-HC" in snapshot.bonds.types

        forces, _ = to_hoomd_forcefield(top=top, r_cut=1.4, base_units=base_units)
        assert "CT-CT-CT-HC" in list(forces["dihedrals"][0].params)
        for conntype in snapshot.dihedrals.types:
            assert conntype in list(forces["dihedrals"][0].params)
