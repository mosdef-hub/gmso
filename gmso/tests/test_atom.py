import numpy as np
import pytest
import unyt as u
from pydantic import ValidationError

from gmso.core.atom import Atom
from gmso.core.atom_type import AtomType
from gmso.core.element import Lithium, Sulfur
from gmso.exceptions import GMSOError
from gmso.tests.base_test import BaseTest


class TestSite(BaseTest):
    def test_new_site(self):
        atom = Atom(name="atom")
        assert atom.name == "atom"
        assert np.all(np.isnan(atom.position))
        assert atom.position.units == u.nm

    def test_dtype(self):
        atom = Atom(name="atom", position=u.nm * np.zeros(3))
        assert atom.position.dtype == float
        assert isinstance(atom.position, u.unyt_array)
        assert isinstance(atom.position, np.ndarray)

    def test_name_none(self):
        atom = Atom()
        assert atom.name == "Atom"

    def test_setters_and_getters(self):
        atom = Atom(name="Lithium", element=Lithium, charge=1, mass=6.941)

        assert atom.name == "Lithium"
        assert atom.element == Lithium
        assert atom.charge == 1 * u.elementary_charge
        assert atom.mass == 6.941 * u.gram / u.mol

        atom.name = "Sulfur"
        atom.element = Sulfur
        atom.charge = -1
        atom.mass = 32.065

        assert atom.name == "Sulfur"
        assert atom.element == Sulfur
        assert atom.charge == -1 * u.elementary_charge
        assert atom.mass == 32.065 * u.gram / u.mol

        atom.mass = None
        atom.charge = None

    def test_set_with_atom_type(self):
        lithium_type = AtomType(mass=6.941, charge=1)
        atom = Atom(name="Lithium")
        atom.atom_type = lithium_type

        assert atom.charge == 1 * u.elementary_charge
        assert atom.mass == 6.941 * u.gram / u.mol

    @pytest.mark.parametrize(
        "position",
        [
            [0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [[0.0, 0.0], [0.0, 0.0]],
            ["a", "b", "c"],
            ["a", 1, 1],
        ],
    )
    def test_bad_pos_input(self, position):
        with pytest.raises((u.exceptions.InvalidUnitOperation, ValueError)):
            Atom(name="atom", position=u.nm * position)

    def test_pos_restraint(self, spce_water):
        for site in spce_water.sites:
            for ax in ["kx", "ky", "kz"]:
                assert site.restraint[ax] == 1000 * u.Unit("kJ/(mol*nm**2)")
                site.restraint[ax] = 500 * u.Unit("kJ/(mol*nm**2)")

    def test_equivalence(self):
        ref = Atom(name="atom", position=u.nm * np.zeros(3))
        same_atom = Atom(name="atom", position=u.nm * np.zeros(3))
        other_pos = Atom(name="atom", position=u.nm * np.ones(3))
        other_name = Atom(name="atom", position=u.nm * np.ones(3))
        # Two sites are never equivalent
        assert ref != same_atom
        assert ref != other_pos
        assert ref != other_name

    def test_atom_tuple_position(self):
        atom = Atom()
        atom.position = (2.0, 2.0, 2.0)
        u.assert_allclose_units(atom.position, [2.0, 2.0, 2.0] * u.nm)

    def test_atom_list_position(self):
        atom = Atom()
        atom.position = [2.0, 2.0, 2.0]
        u.assert_allclose_units(atom.position, (2.0, 2.0, 2.0) * u.nm)

    def test_atom_invalid_position(self):
        atom = Atom(name="invalidSite")
        with pytest.raises(ValueError) as e:
            atom.position = [0, 0, 0, 0]
            assert (
                "Position of shape (3,) is not valid. "
                "Accepted values: (a.) 3-tuple, (b.) list of length 3 "
                "(c.) np.array or unyt.unyt_array of shape (3,)" in e
            )

    def test_position_assignment_invalid(self):
        atom1 = Atom(name="Site")

        with pytest.raises(ValidationError) as e:
            atom1.position = "invalid"
