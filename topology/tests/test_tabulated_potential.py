import numpy as np
import unyt as u
import pytest

from topology.core.tabulated_potential import TabulatedPotential
from topology.tests.base_test import BaseTest
from topology.utils.testing import allclose
from topology.exceptions import TopologyError


class TestPotential(BaseTest):
    def test_new_potential(self):
        pot = TabulatedPotential(
            r = np.linspace(0, 1),
            v = np.linspace(5, 0),
        )

    def test_getters_setters(self, tab_pot):
        assert np.allclose(tab_pot.r, np.linspace(0, 1))
        assert np.allclose(tab_pot.v, np.linspace(5, 0))

        tab_pot.r = np.linspace(0, 10)
        tab_pot.v = np.linspace(2, 0)

        assert np.allclose(tab_pot.r, np.linspace(0, 10))
        assert np.allclose(tab_pot.v, np.linspace(2, 0))

    def test_set_both(self, tab_pot):
        tab_pot.set_r_and_v(
            r=np.linspace(0, 1, 40),
            v=np.linspace(5, 0, 40),
        )

        assert np.allclose(tab_pot.r, np.linspace(0, 1, 40))
        assert np.allclose(tab_pot.v, np.linspace(5, 0, 40))

    def test_bad_setters(self, tab_pot):
        with pytest.raises(TopologyError):
            tab_pot.r = np.linspace(0, 10, num=73)
        with pytest.raises(TopologyError):
            tab_pot.v = np.eye(4)
