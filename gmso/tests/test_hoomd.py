import forcefield_utilities as ffutils
import hoomd
import numpy as np
import pytest
import unyt as u
from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield

from gmso.external import from_mbuild
from gmso.external.convert_hoomd import to_hoomd_forcefield, to_hoomd_snapshot
from gmso.parameterization import apply
from gmso.tests.base_test import BaseTest
from gmso.utils.io import has_hoomd, has_mbuild, import_

if has_hoomd:
    hoomd = import_("hoomd")
if has_mbuild:
    mb = import_("mbuild")


@pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
@pytest.mark.skipif(not has_mbuild, reason="mbuild not installed")
class TestGsd(BaseTest):
    def test_mbuild_comparison(self):
        compound = mb.load("CCC", smiles=True)
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=20)
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

        mb_snapshot, mb_forcefield, ref_vals = create_hoomd_forcefield(
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
            if not isinstance(mb_force, hoomd.md.long_range.pppm.Coulomb):
                keys = mb_force.params.param_dict.keys()
                for key in keys:
                    mb_params = mb_force.params.param_dict[key]
                    gmso_params = gmso_force.params.param_dict[key]
                    variables = mb_params.keys()
                    for var in variables:
                        assert np.isclose(mb_params[var], gmso_params[var])

    def test_hoomd_simulation(self):
        compound = mb.load("CCC", smiles=True)
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=200)
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

    def test_hoomd_simulation_auto_scaled(self):
        compound = mb.load("CCC", smiles=True)
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=200)
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
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=100)
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
        com_box = mb.packing.fill_box(compound, box=[5, 5, 5], n_compounds=100)
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
            base_units=base_units,
            pppm_kwargs={"resolution": (64, 64, 64), "order": 7},
        )
