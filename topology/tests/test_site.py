import numpy as np
import unyt as u
import pytest

from topology.core.site import Site
from topology.tests.base_test import BaseTest


class TestSite(BaseTest):
    def test_new_site(self):
        site = Site(name='site')
        assert site.name == 'site'

    def test_dtype(self):
        site = Site(name='site', position=u.nm*np.zeros(3))
        assert site.position.dtype == float
        assert isinstance(site.position, u.unyt_array)
        assert isinstance(site.position, np.ndarray)

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

        assert ref == same_site
        assert ref != other_pos
        assert ref != other_name
