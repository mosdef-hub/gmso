import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

from gmso import Topology
from gmso.core.atom import Atom
from gmso.core.box import Box
from gmso.external.convert_parmed import from_parmed
from gmso.formats.gro import _prepare_atoms, read_gro, write_gro
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path
from gmso.utils.io import get_fn, has_mbuild, has_parmed, import_

if has_parmed:
    pmd = import_("parmed")

if has_mbuild:
    mb = import_("mbuild")


@pytest.mark.skipif(not has_parmed, reason="ParmEd is not installed")
class TestGro(BaseTest):
    def test_read_gro(self):
        top = Topology.load(get_fn("acn.gro"))

        assert top.name == "ACN"
        assert top.n_sites == 6
        assert_allclose_units(
            top.box.lengths, 4 * np.ones(3) * u.nm, rtol=1e-5, atol=1e-8
        )

        top = Topology.load(get_fn("350-waters.gro"))

        assert top.name == "Generic title"
        assert top.n_sites == 1050
        assert_allclose_units(
            top.box.lengths, 2.20866 * np.ones(3) * u.nm, rtol=1e-5, atol=1e-8
        )

    def test_wrong_n_atoms(self):
        with pytest.raises(ValueError):
            Topology.load(get_fn("too_few_atoms.gro"))
        with pytest.raises(ValueError):
            Topology.load(get_fn("too_many_atoms.gro"))

    def test_write_gro(self):
        top = from_parmed(pmd.load_file(get_fn("ethane.gro"), structure=True))
        top.save("out.gro")

    def test_write_gro_with_shift_coord(self):
        top = from_parmed(pmd.load_file(get_fn("ethane.mol2"), structure=True))
        top.save("out.gro", shift_coord=True)

        read_top = Topology.load("out.gro")
        assert np.all(list(map(lambda x: x.position >= 0, read_top.sites)))

    def test_write_gro_non_orthogonal(self):
        top = from_parmed(pmd.load_file(get_fn("ethane.gro"), structure=True))
        top.box.angles = u.degree * [90, 90, 120]
        top.save("out.gro")

    @pytest.mark.skipif(not has_mbuild, reason="mBuild not installed.")
    def test_benzene_gro(self):
        import mbuild as mb
        from mbuild.packing import fill_box

        from gmso.external import from_mbuild

        benzene = mb.load(get_fn("benzene.mol2"))
        benzene.children[0].name = "Benzene"
        box_of_benzene = fill_box(compound=benzene, n_compounds=5, density=1)
        top = from_mbuild(box_of_benzene)
        top.save("benzene.gro")

        reread = Topology.load("benzene.gro")
        for site, ref_site in zip(reread.sites, top.sites):
            assert site.molecule.name == ref_site.molecule.name[:5]
            assert site.molecule.number == ref_site.molecule.number

    def test_prepare_atoms(self):
        top = Topology()
        ref = Atom(name="atom1", position=[0.0, 0.0, 3.0], molecule=("mol", 0))
        top.add_site(ref)

        line = _prepare_atoms(top, top.positions, 5)
        assert line == "    1mol  atom1    1   0.00000   0.00000   3.00000\n"

        top = Topology()
        ref = Atom(name="atom", position=[0.0, 0.0, 0.0])
        top.add_site(ref)

        line = _prepare_atoms(top, top.positions, 5)
        assert line == "    1MOL   atom    1   0.00000   0.00000   0.00000\n"

    @pytest.mark.skipif(not has_mbuild, reason="mBuild not installed.")
    def test_resid_for_mol(self):
        # test adding different molecules to the system

        import mbuild as mb

        from gmso.external import from_mbuild

        ethane = mb.lib.molecules.Ethane()
        methane = mb.lib.molecules.Methane()
        system = mb.Compound()

        system.add(mb.clone(ethane))
        system.add(mb.clone(ethane))
        system.add(mb.clone(methane))
        system.add(mb.clone(methane))

        top = from_mbuild(system)
        top.save("ethane_methane.gro")

        reread = Topology.load("ethane_methane.gro")
        nums = set([site.molecule.number for site in reread.sites])
        assert nums == {0, 1, 2, 3}

    def test_no_mol_name(self):
        # here we will just add sites with no molecule information to
        # ensure that residues are labeled with the same mol number

        top = Topology()
        for i in range(0, 2):
            ref = Atom(name="atom", position=[0.0, 0.0, 0.0])
            top.add_site(ref)
        box = Box(2 * u.nm * np.ones(3))
        top.box = box
        top.save("temp_system.gro")
        reread = Topology.load("temp_system.gro")
        nums = set([site.molecule.number for site in reread.sites])
        assert nums == {0}

    def test_res_naming(self):
        top = Topology()
        ref = Atom(
            name="mol_atom", position=[0.0, 0.0, 0.0], molecule=("test", 0)
        )
        top.add_site(ref)

        for i in range(0, 2):
            ref = Atom(name="atom", position=[0.0, 0.0, 0.0])
        top.add_site(ref)
        box = Box(2 * u.nm * np.ones(3))

        top.box = box
        top.save("temp1.gro", overwrite=True)

        reread = Topology.load("temp1.gro")
        nums = set([site.molecule.number for site in reread.sites])
        assert nums == {0, 1}

        top = Topology()
        ref = Atom(
            name="mol_atom", position=[0.0, 0.0, 0.0], molecule=("test", 0)
        )
        top.add_site(ref)
        ref = Atom(
            name="mol_atom", position=[0.0, 0.0, 0.0], molecule=("test", 0)
        )
        top.add_site(ref)

        ref = Atom(
            name="mol_atom", position=[0.0, 0.0, 0.0], molecule=("test", 1)
        )
        top.add_site(ref)
        ref = Atom(
            name="mol_atom", position=[0.0, 0.0, 0.0], molecule=("test", 1)
        )
        top.add_site(ref)

        for i in range(0, 2):
            ref = Atom(name="atom", position=[0.0, 0.0, 0.0])
            top.add_site(ref)
        box = Box(2 * u.nm * np.ones(3))
        top.box = box
        top.save("temp2.gro", overwrite=True)

        reread = Topology.load("temp2.gro")
        nums = set([site.molecule.number for site in reread.sites])
        assert nums == {0, 1, 2}

        top = Topology()
        ref = Atom(
            name="mol_atom", position=[0.0, 0.0, 0.0], molecule=("test", 0)
        )
        top.add_site(ref)
        ref = Atom(
            name="mol_atom", position=[0.0, 0.0, 0.0], molecule=("test", 0)
        )
        top.add_site(ref)

        ref = Atom(name="resA", position=[0.0, 0.0, 0.0], residue=("resA", 0))
        top.add_site(ref)
        ref = Atom(name="resB", position=[0.0, 0.0, 0.0], residue=("resB", 1))
        top.add_site(ref)

        for i in range(0, 2):
            ref = Atom(name="atom", position=[0.0, 0.0, 0.0])
            top.add_site(ref)
        box = Box(2 * u.nm * np.ones(3))
        top.box = box
        top.save("temp3.gro", overwrite=True)

        reread = Topology.load("temp3.gro")
        nums = set([site.molecule.number for site in reread.sites])
        assert nums == {0, 1, 2, 3}

    @pytest.mark.parametrize("fixture", ["benzene_ua_box", "benzene_aa_box"])
    def test_full_loop_gro_molecule(self, fixture, request):
        top = request.getfixturevalue(fixture)
        top.save("benzene.gro")

        # Re-read in and compare with reference
        top = Topology.load("benzene.gro")

        refs = {
            "benzene_aa_box": "benzene.gro",
            "benzene_ua_box": "restrained_benzene_ua.gro",
        }
        ref = Topology.load(get_path(refs[fixture]))

        assert len(top.sites) == len(ref.sites)
        assert top.unique_site_labels("molecule") == ref.unique_site_labels(
            "molecule"
        )

        if top == "benzene_ua_box":
            for site in top.sites:
                assert site.molecule.name == "Com"
                assert site.name == "_CH"
        elif top == "benzene_aa_box":
            for site in top.sites:
                assert site.molecule.name == "Benze"
                assert site.name in ["C", "H"]
