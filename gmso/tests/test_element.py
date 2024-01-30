import numpy as np
import pytest
import unyt as u

from gmso.core import element
from gmso.core.atom_type import AtomType
from gmso.core.element import Carbon
from gmso.exceptions import GMSOError
from gmso.tests.base_test import BaseTest


class TestElement(BaseTest):
    def test_element(self):
        carbon = Carbon

        assert carbon.name == "carbon"
        assert carbon.symbol == "C"
        assert carbon.mass == element.Carbon.mass

    def test_element_by_name(self):
        for idx, name in enumerate(["Carbon", "carbon", " CarBon 12 "]):
            if idx == 1:
                carbon = element.element_by_name(name, verbose=True)
            else:
                with pytest.warns(UserWarning):
                    carbon = element.element_by_name(name, verbose=True)
            assert carbon.name == element.Carbon.name
            assert carbon.symbol == element.Carbon.symbol
            assert carbon.mass == element.Carbon.mass

    def test_element_by_symbol(self):
        for idx, symbol in enumerate(["N", "n", " N7"]):
            if idx == 0:
                nitrogen = element.element_by_symbol(symbol, verbose=True)
            else:
                with pytest.warns(UserWarning):
                    nitrogen = element.element_by_symbol(symbol, verbose=True)
            assert nitrogen.name == element.Nitrogen.name
            assert nitrogen.symbol == element.Nitrogen.symbol
            assert nitrogen.mass == element.Nitrogen.mass

    def test_element_by_atomic_number(self):
        for number in [8, "8", "08", "Oxygen-08"]:
            oxygen = element.element_by_atomic_number(number)

            assert oxygen.name == element.Oxygen.name
            assert oxygen.symbol == element.Oxygen.symbol
            assert oxygen.mass == element.Oxygen.mass

    def test_element_by_mass(self):
        for mass in ["Fluorine-19", 19, 19 * u.amu]:
            fluorine = element.element_by_mass(mass)

            assert fluorine.name == element.Fluorine.name
            assert fluorine.symbol == element.Fluorine.symbol
            assert fluorine.mass == element.Fluorine.mass

    def test_element_by_mass_exactness(self):
        cobalt = element.element_by_mass(58.9)
        nickel = element.element_by_mass(58.7)
        chlorine = element.element_by_mass(35, exact=False)

        assert cobalt == element.Cobalt
        assert nickel == element.Nickel
        assert chlorine == element.Chlorine

    def test_element_by_smarts_string(self):
        SMARTS = [
            "C",
            "[C;X3](C)(O)C",
            "[C;X4](C)(C)(C)H",
            "[C;X4;r3]1(H)(H)[C;X4;r3][C;X4;r3]1",
            "[C;X3]([O;X1])[N;X3]",
            "[C;X3;r6]1(OH)[C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6]1",
        ]

        for smarts in SMARTS:
            carbon = element.element_by_smarts_string(smarts)

            assert carbon.name == element.Carbon.name
            assert carbon.symbol == element.Carbon.symbol
            assert carbon.mass == element.Carbon.mass

    def test_element_by_atom_type(self):
        carbon_type = AtomType(mass=12.011, definition="C", name="C")
        mass_only = AtomType(mass=12.011)
        def_only = AtomType(definition="C")
        name_only = AtomType(name="C")

        for atom_type in [carbon_type, mass_only, def_only, name_only]:
            carbon = element.element_by_atom_type(atom_type)

            assert carbon.name == element.Carbon.name
            assert carbon.symbol == element.Carbon.symbol
            assert carbon.mass == element.Carbon.mass

    def test_all_elements(self):
        for num in range(1, 119):
            elem = element.element_by_atomic_number(num)
            assert elem.atomic_number == num
            assert element.__hash__

    @pytest.mark.parametrize("atomic_number", [0, -1, 1000, 17.02])
    def test_bad_atomic_number(self, atomic_number):
        with pytest.raises(GMSOError):
            element.element_by_atomic_number(atomic_number)
