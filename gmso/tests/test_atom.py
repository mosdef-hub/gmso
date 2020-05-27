import numpy as np
import unyt as u
import pytest

from pydantic import ValidationError

from gmso.core.atom import Atom
from gmso.core.bond import Bond
from gmso.core.element import Lithium, Sulfur
from gmso.core.atom_type import AtomType
from gmso.tests.base_test import BaseTest
from gmso.exceptions import GMSOError
from gmso.utils.testing import allclose


class TestAtom(BaseTest):
    def test_new_atom(self):
        atom = Atom(name='atom')
        assert atom.name == 'atom'
        assert atom.position is None

    def test_dtype(self):
        atom = Atom(name='atom', position=u.nm * np.zeros(3))
        assert atom.position.dtype == float
        assert isinstance(atom.position, u.unyt_array)
        assert isinstance(atom.position, np.ndarray)

    def test_name_none(self):
        atom = Atom(name=None)
        assert atom.name == 'Atom'

    def test_setters_and_getters(self):
        atom = Atom(name='Lithium', element=Lithium, charge=1, mass=6.941)

        assert atom.name == 'Lithium'
        assert atom.element == Lithium
        assert atom.charge == 1*u.elementary_charge
        assert atom.mass == 6.941*u.gram/u.mol

        atom.name = 'Sulfur'
        atom.element = Sulfur
        atom.charge = -1
        atom.mass=32.065

        assert atom.name == 'Sulfur'
        assert atom.element == Sulfur
        assert atom.charge == -1*u.elementary_charge
        assert atom.mass == 32.065*u.gram/u.mol

    def test_set_with_atom_type(self):
        lithium_type = AtomType(mass=6.941, charge=1)
        atom = Atom(name='Lithium')
        atom.atom_type = lithium_type

        assert atom.charge == 1*u.elementary_charge
        assert atom.mass == 6.941*u.gram/u.mol

    @pytest.mark.parametrize('position', [[0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0], [[0.0, 0.0], [0.0, 0.0]],
        ['a', 'b', 'c'], ['a', 1, 1]])
    def test_bad_pos_input(self, position):
        with pytest.raises((u.exceptions.InvalidUnitOperation, ValueError)):
            Atom(name='atom', position=u.nm * position)

    def test_equivalence(self):
        ref = Atom(name='atom', position=u.nm * np.zeros(3))
        same_atom = Atom(name='atom', position=u.nm * np.zeros(3))
        other_pos = Atom(name='atom', position=u.nm * np.ones(3))
        other_name = Atom(name='atom', position=u.nm * np.ones(3))
        # Two atoms are never equivalent
        assert ref != same_atom
        assert ref != other_pos
        assert ref != other_name

    def test_position_assignment(self):
        atom1 = Atom(name='Atom', position=[1.0, 1.0, 1.0])
        new_position = np.array([2., 2., 2.])
        atom1.position = new_position
        assert allclose(atom1.position, u.unyt_array(new_position, u.nm))

    def test_position_assignment_invalid(self):
        atom1 = Atom(name='Atom')
        with pytest.raises(ValidationError) as e:
            atom1.position = 'invalid'
            assert "Converting object of type <class 'str'> failed with error Tried to multiply a Unit " \
                   "object with 'a' (type <class 'str'>). This behavior is undefined." in e
