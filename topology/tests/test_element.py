import numpy as np

from topology.core.element import Carbon
from topology.tests.base_test import BaseTest


class TestElement(BaseTest):
    def test_element(self):
        carbon = Carbon

        assert carbon.name == 'carbon'
        assert carbon.symbol == 'C'
        assert np.isclose(carbon.mass, 12.011)
