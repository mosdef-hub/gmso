import pytest

from gmso.abc.abstract_site import Molecule, Residue, Site
from gmso.tests.base_test import BaseTest


class TestSite(BaseTest):
    def test_molecule(self):
        molecule = Molecule()
        molecule.__repr__()
        assert molecule.number == 0
        assert molecule.isrigid is False

    def test_residue(self):
        residue = Residue()
        residue.__repr__()
        assert residue.number == 0

    def test_site(self):
        with pytest.raises(TypeError):
            Site()
