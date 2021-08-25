import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

from gmso import Topology
from gmso.formats.mol2 import from_mol2
from gmso.tests.base_test import BaseTest
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

        match = "No charges were detected for site C with index 1"
        with pytest.warns(UserWarning) as record:
            top = Topology.load(get_fn("ethane.mol2"))
        assert record[3].message.args[0] == match
        assert list(top.sites)[0].charge is None

        match = "No element detected for site C with index1, consider manually adding the element to the topology"
        with pytest.warns(UserWarning) as record:
            Topology.load(get_fn("benzene.mol2"))
        assert record[3].message.args[0] == match

    def test_residue(self):
        top = Topology.load(get_fn("ethanol_aa.mol2"))
        assert np.all([site.residue_label == 'ETO' for site in top.sites])
        assert np.all([site.residue_index == 1 for site in top.sites])

        top = Topology.load(get_fn("benzene_ua.mol2"), site_type="lj")
        assert np.all([site.residue_label == 'BEN1' for site in top.iter_sites("residue_label","BEN1")])
        assert np.all([site.residue_index == 1 for site in top.iter_sites("residue_label","BEN1")])
        assert np.all([site.residue_label == 'BEN2' for site in top.iter_sites("residue_label","BEN2")])
        assert np.all([site.residue_index == 2 for site in top.iter_sites("residue_label","BEN2")])

    def test_lj_system(self):
        top = Topology.load(get_fn("methane.mol2"), site_type="lj")
        assert np.all([site.element == None for site in top.sites])

    def test_wrong_path(self):
        with pytest.raises(OSError) as e:
            Topology.load("not_a_file.mol2")
        assert e.value.args[0] == "Provided path to file that does not exist"
        top = Topology.load(get_fn("ethanegro.mol2"))
        assert len(top.sites) == 0
        assert len(top.bonds) == 0

    def test_broken_files(self):
        match = "The record type indicator @<TRIPOS>MOLECULE_extra_text\n is not supported"
        with pytest.warns(UserWarning) as record:
            Topology.load(get_fn("broken.mol2"))
        assert record[0].message.args[0] == match
        match = "This mol2 file has two boxes to be read in, only reading in one with dimensions Box(a=0.72 nm, b=0.9337999999999999 nm, c=0.6646 nm, alpha=90.0 degree, beta=90.0 degree, gamma=90.0 degree)"
        assert record[90].message.args[0] == match
