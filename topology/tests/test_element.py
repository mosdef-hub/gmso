import pytest
import numpy as np
import unyt as u

from topology.core import element
from topology.core.element import Carbon
from topology.tests.base_test import BaseTest
from topology.exceptions import TopologyError


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
        for number in [8, '8', '08', 'Oxygen-08']:
            oxygen = element.element_by_atomic_number(number)

            assert oxygen.name == element.Oxygen.name
            assert oxygen.symbol == element.Oxygen.symbol
            assert oxygen.mass == element.Oxygen.mass

        #Test element_by_mass
        for mass in ['Fluorine-19', 19, 19 * u.amu]:
            fluorine = element.element_by_mass(mass)

            assert fluorine.name == element.Fluorine.name
            assert fluorine.symbol == element.Fluorine.symbol
            assert fluorine.mass == element.Fluorine.mass 

        for smarts in ['C', '[C;X3](C)(O)C', '[C;X4](C)(C)(C)H']:
            carbon = element.element_by_smarts_string(smarts)

            assert carbon.name == element.Carbon.name
            assert carbon.symbol == element.Carbon.symbol
            assert carbon.mass == element.Carbon.mass

        #Additional element_by_mass test
        cobalt = element.element_by_mass(58.9)
        nickel = element.element_by_mass(58.7)
        chlorine = element.element_by_mass(35, exact=False)

        assert cobalt == element.Cobalt
        assert nickel == element.Nickel
        assert chlorine == element.Chlorine

    def test_all_elements(self):
        for num in range(1, 119):
            elem = element.element_by_atomic_number(num)
            assert elem.atomic_number == num

    @pytest.mark.parametrize('atomic_number', [0, -1, 1000, 17.02])
    def test_bad_atomic_number(self, atomic_number):
        with pytest.raises(TopologyError):
            element.element_by_atomic_number(atomic_number)
