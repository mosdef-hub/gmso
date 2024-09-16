from gmso import ForceField, Topology
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
        ptop = apply(
            top, ff, identify_connections=True, ignore_params=["dihedral", "improper"]
        )
        ptop.save("charmm_dihedral.prm")

    def test_par_parameters(self, ethane, oplsaa_forcefield):
        from gmso.parameterization import apply

        ptop = apply(ethane, oplsaa_forcefield)

        ptop.save(filename="ethane-opls.par")
        from parmed.charmm import CharmmParameterSet

        pset = CharmmParameterSet.load_set(pfile="ethane-opls.par")
        assert len(pset.bond_types) == 3
        assert len(pset.angle_types) == 3
        assert len(pset.atom_types) == 2

    def test_parameter_forms(self, ethane):
        # test that the values are correct
        # write mbuild file and compare to gmso file
        pass

    def test_raise_errors(self):
        # only takes harmonic bonds
        # only takes harmonic of ub angles
        pass
