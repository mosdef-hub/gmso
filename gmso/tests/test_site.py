import numpy as np
import unyt as u
import pytest

from gmso.core.site import Site
from gmso.core.bond import Bond
from gmso.core.element import Lithium, Sulfur
from gmso.core.atom_type import AtomType
from gmso.tests.base_test import BaseTest
from gmso.exceptions import GMSOError


class TestSite(BaseTest):
    def test_new_site(self):
        site = Site(name='site')
        assert site.name == 'site'

    def test_dtype(self):
        site = Site(name='site', position=u.nm*np.zeros(3))
        assert site.position.dtype == float
        assert isinstance(site.position, u.unyt_array)
        assert isinstance(site.position, np.ndarray)

    def test_name_none(self):
        site = Site(name=None)
        assert site.name == 'Site'

    def test_setters_and_getters(self):
        site = Site(name='Lithium', element=Lithium, charge=1, mass=6.941)

        assert site.name == 'Lithium'
        assert site.element == Lithium
        assert site.charge == 1*u.elementary_charge
        assert site.mass == 6.941*u.gram/u.mol

        site.name = 'Sulfur'
        site.element = Sulfur
        site.charge = -1
        site.mass=32.065

        assert site.name == 'Sulfur'
        assert site.element == Sulfur
        assert site.charge == -1*u.elementary_charge
        assert site.mass == 32.065*u.gram/u.mol

        site.element, site.charge, site.mass = None, None, None
        assert site.element is None
        assert site.charge is None
        assert site.mass is None

    def test_set_with_atom_type(self):
        lithium_type = AtomType(mass=6.941, charge=1)
        site = Site(name='Lithium')
        site.atom_type = lithium_type

        assert site.charge == 1*u.elementary_charge
        assert site.mass == 6.941*u.gram/u.mol

    @pytest.mark.parametrize('position', [[0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0], [[0.0, 0.0], [0.0, 0.0]],
        ['a', 'b', 'c'], ['a', 1, 1]])
    def test_bad_pos_input(self, position):
        with pytest.raises((u.exceptions.InvalidUnitOperation, ValueError)):
            Site(name='site', position=u.nm*position)

    def test_equivalence(self):
        ref = Site(name='site', position=u.nm*np.zeros(3))
        same_site = Site(name='site', position=u.nm*np.zeros(3))
        other_pos = Site(name='site', position=u.nm*np.ones(3))
        other_name = Site(name='site', position=u.nm*np.ones(3))
        # Two sites are never equivalent
        assert ref != same_site
        assert ref != other_pos
        assert ref != other_name
