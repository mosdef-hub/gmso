import pytest

from gmso import ForceField, Topology
from gmso.exceptions import EngineIncompatibilityError
from gmso.parameterization import apply
from gmso.tests.base_test import BaseTest
from gmso.utils.io import get_fn

# TODO: Sorting of atomtypes info in connection types
# TODO: Test of correct potential forms
# TODO: Handle iterating over unique types
# TODO: Handle multiple values in charmm dihedral k as np array


class TestPar(BaseTest):
    def test_save_charmm(self):
        top = Topology.load(get_fn("charmm_dihedral.mol2"))
        for site in top.sites:
            site.name = f"_{site.name}"

        ff = ForceField(get_fn("charmm_truncated.xml"))
        ptop = apply(top, ff, identify_connections=True, ignore_params=["improper"])
        ptop.save("charmm_dihedral.prm")

    def test_par_parameters(self):
        from gmso.parameterization import apply

        top = Topology.load(get_fn("charmm_improper.mol2"))
        ff = ForceField(get_fn("charmm_amoeba.xml"))
        for site in top.sites:
            site.name = f"_{site.name}"

        ptop = apply(
            top, ff, identify_connections=True, ignore_params=["improper", "dihedral"]
        )

        ptop.save(filename="ethane-opls.par")
        from parmed.charmm import CharmmParameterSet

        pset = CharmmParameterSet.load_set(pfile="ethane-opls.par")
        assert len(pset.bond_types) == 10  # each bond counted twice
        assert len(pset.angle_types) == 10  # each angle counted twice
        assert len(pset.dihedral_types) == 2
        assert len(pset.improper_types) == 1
        assert len(pset.atom_types) == 4

    def test_prm_incompatibile_types(self, ethane, oplsaa_forcefield):
        from gmso.parameterization import apply

        ptop = apply(ethane, oplsaa_forcefield, identify_connections=True)
        with pytest.raises(EngineIncompatibilityError):
            ptop.save(filename="ethane-opls.par")

    def test_parameter_forms(self, ethane):
        # test that the values are correct
        # write mbuild file and compare to gmso file
        pass

    def test_raise_errors(self):
        # only takes harmonic bonds
        # only takes harmonic of ub angles
        pass
