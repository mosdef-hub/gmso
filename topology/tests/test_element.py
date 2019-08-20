import numpy as np

from topology.core import element
from topology.core.element import Carbon
from topology.tests.base_test import BaseTest
import unyt as u

class TestElement(BaseTest):
    def test_element(self):
        carbon = Carbon

        assert carbon.name == 'carbon'
        assert carbon.symbol == 'C'
        assert carbon.mass == element.Carbon.mass

    def test_element_by(self):
        #Test element_by_name
        for name in ['Carbon', 'carbon', ' CarBon 12 ']:
            carbon = element.element_by_name(name)

            assert carbon.name == element.Carbon.name
            assert carbon.symbol == element.Carbon.symbol
            assert carbon.mass == element.Carbon.mass

        #Test element_by_symbol
        for symbol in ['N', 'n', ' N7']:
            nitrogen = element.element_by_symbol(symbol)

            assert nitrogen.name == element.Nitrogen.name
            assert nitrogen.symbol == element.Nitrogen.symbol
            assert nitrogen.mass ==  element.Nitrogen.mass

        #Test element_by_atomic_number
        for number in [8, '8', '08', 'Oxygen08']:
            oxygen = element.element_by_atomic_number(number)

            assert oxygen.name == element.Oxygen.name
            assert oxygen.symbol == element.Oxygen.symbol
            assert oxygen.mass == element.Oxygen.mass
