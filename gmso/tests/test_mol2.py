import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

from gmso import Topology
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path
from gmso.utils.io import get_fn


class TestMol2(BaseTest):
    def test_read_mol2(self):
        top = Topology.load(get_fn("parmed.mol2"))
        assert top.name == "parmed"
        assert top.n_sites == 8
        assert_allclose_units(
            top.box.lengths,
            ([8.2693, 7.9100, 6.6460] * u.Å).to("nm"),
            rtol=1e-5,
            atol=1e-8,
        )
        assert list(top.sites)[0].element.name == "carbon"
        assert_allclose_units(
            list(top.sites)[0].element.mass,
            np.array(1.9944733e-26) * u.kg,
            rtol=1e-5,
            atol=1e-8,
        )

        top = Topology.load(get_fn("tip3p.mol2"))
        assert top.name == "tip3p"
        assert top.n_sites == 3
        assert_allclose_units(
            top.box.lengths, 3.0130 * np.ones(3) * u.Å, rtol=1e-5, atol=1e-8
        )
        positions_check = [
            [0.061, 0.1, 0.1],
            [0.017, 0.09, 0.177],
            [0.011, 0.154, 0.04],
        ]
        for check, site in zip(positions_check, top.sites):
            assert_allclose_units(
                site.position,
                check * u.nm,
                rtol=1e-5,
                atol=1e-8,
            )

        top = Topology.load(get_fn("vmd.mol2"))
        assert top.name == "vmd"
        assert top.n_sites == 6
        assert len(top.bonds) == 5
        assert top.bonds[0].connection_members[0] == top.sites[0]
        assert top.box == None

        with pytest.warns(
            UserWarning,
            match=r"No charge was detected for site C with index 1",
        ):
            top = Topology.load(get_fn("ethane.mol2"), verbose=True)
        assert list(top.sites)[0].charge is None

    def test_residue(self):
        top = Topology.load(get_fn("ethanol_aa.mol2"))
        assert np.all([site.residue[0] == "ETO" for site in top.sites])
        assert np.all([site.residue[1] == 1 for site in top.sites])

        top = Topology.load(get_fn("benzene_ua.mol2"), site_type="lj")
        assert np.all(
            [
                site.residue[0] == "BEN1"
                for site in top.iter_sites_by_residue("BEN1")
            ]
        )
        assert np.all(
            [site.residue[1] == 1 for site in top.iter_sites_by_residue("BEN1")]
        )
        assert np.all(
            [
                site.residue[0] == "BEN2"
                for site in top.iter_sites_by_residue("BEN2")
            ]
        )
        assert np.all(
            [site.residue[1] == 2 for site in top.iter_sites_by_residue("BEN2")]
        )

    def test_lj_system(self):
        top = Topology.load(get_fn("methane.mol2"), site_type="lj")
        assert np.all([site.element == None for site in top.sites])

    def test_no_charge_lj(self):
        with pytest.warns(
            UserWarning,
            match=r"No charge was detected for site .* with index \d+$",
        ):
            top = Topology.load(
                get_path("methane_missing_charge.mol2"),
                site_type="lj",
                verbose=True,
            )

    def test_wrong_path(self):
        with pytest.raises(
            OSError, match=r"Provided path to file that does not exist"
        ):
            Topology.load("not_a_file.mol2")
        top = Topology.load(get_fn("ethanegro.mol2"))
        assert len(top.sites) == 0
        assert len(top.bonds) == 0

    def test_broken_files(self):
        with pytest.warns(
            UserWarning,
            match=r"The record type indicator @<TRIPOS>MOLECULE_extra_text is not supported. Skipping current section and moving to the next RTI header.",
        ):
            Topology.load(get_fn("broken.mol2"))
        with pytest.warns(
            UserWarning,
            match=r"This mol2 file has two boxes to be read in, only reading in one with dimensions Box\(a=0.72",
        ):
            Topology.load(get_fn("broken.mol2"), verbose=True)

    def test_benzene_mol2_elements(self):
        top = Topology.load(get_fn("benzene.mol2"))

        for atom in top.sites:
            assert atom.element.name in {"hydrogen", "carbon"}

    def test_neopentane_mol2_elements(self):
        with pytest.warns(
            UserWarning,
            match=r"No element detected for site .+ with index \d+, "
            r"consider manually adding the element to the topology$",
        ):
            top = Topology.load(get_fn("neopentane.mol2"), verbose=True)

    def test_mol2_residues(self):
        top = Topology.load(get_fn("parmed.mol2"))
        assert np.all(
            np.array([site.residue.name for site in top.sites]) == "RES"
        )
        assert np.all(
            np.array([site.residue.number for site in top.sites]) == 1
        )

    def test_mol2_molecules(self):
        top = Topology.load(get_fn("methane.mol2"))
        assert np.all(
            np.array([site.molecule.name for site in top.sites]) == "MET"
        )
        assert np.all(
            np.array([site.molecule.number for site in top.sites]) == 0
        )

    def test_mol2_group(self):
        # Is there a place to read from mol2 file?
        top = Topology.load(get_fn("ethane.mol2"))
        for site in top.sites:
            site.group = "ethane"
        assert np.all(np.array([site.group for site in top.sites]) == "ethane")
