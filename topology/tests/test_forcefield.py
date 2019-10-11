import pytest

from topology.core.forcefield import Forcefield
from topology.core.atom_type import AtomType
from topology.core.bond_type import BondType
from topology.core.angle_type import AngleType


from topology.exceptions import TopologyError
from topology.tests.base_test import BaseTest


class TestForcefield(BaseTest):
    def test_init(self):
        ff = Forcefield(name='TestFF')
        assert ff.name == 'TestFF'
        assert len(ff.atom_types) == 0
        assert len(ff.bond_types) == 0
        assert len(ff.angle_types) == 0

    def test_setters_getters(self):
        ff = Forcefield(name='TestFF')
        myatomtype = AtomType()
        mybondtype = BondType()
        myangletype = AngleType()

        ff['A'] = myatomtype
        ff['A-B'] = mybondtype
        ff['A-B-C'] = myangletype

        assert len(ff.atom_types) == 1
        assert len(ff.bond_types) == 1
        assert len(ff.angle_types) == 1

        assert ff['A'] == myatomtype
        assert ff['A-B'] == mybondtype
        assert ff['A-B-C'] == myangletype

        assert ff.atom_types['A'] == myatomtype
        assert ff.bond_types['A-B'] == mybondtype
        assert ff.angle_types['A-B-C'] == myangletype

    def test_bad_setters(self):
        ff = Forcefield(name='TestFF')
        myatomtype = AtomType()
        mybondtype = BondType()
        myangletype = AngleType()

        with pytest.raises(TopologyError):
            ff['A'] = mybondtype
        with pytest.raises(TopologyError):
            ff['A-B'] = myangletype
        with pytest.raises(TopologyError):
            ff['A-B-C'] = myatomtype
        with pytest.raises(TopologyError):
            ff['A-B-C-D-E-F'] = 42

    def test_bad_getters(self):
        ff = Forcefield(name='TestFF')
        myatomtype = AtomType()
        mybondtype = BondType()
        myangletype = AngleType()

        ff['A'] = myatomtype
        ff['A-B'] = mybondtype
        ff['A-B-C'] = myangletype

        assert ff['B'] == None
        assert ff['B-C'] == None
        assert ff['B-C-D'] == None
        with pytest.raises(TopologyError):
            ff['A-B-C-D-E-F']

    def test_getters_diff_order(self):
        ff = Forcefield(name='TestFF')
        myatomtype = AtomType()
        mybondtype = BondType()
        myangletype = AngleType()

        ff['A'] = myatomtype
        ff['A-B'] = mybondtype
        ff['A-B-C'] = myangletype

        assert ff['B-A'] == mybondtype
        assert ff['C-B-A'] == myangletype

    def test_getters_wildcards(self):
        ff = Forcefield(name='TestFF')
        atomtype1 = AtomType(name='atomtype1')
        atomtype2 = AtomType(name='atomtype2')
        bondtype1 = BondType(name='bondtype1')
        bondtype2 = BondType(name='bondtype2')
        angletype1 = AngleType(name='angletype1')
        angletype2 = AngleType(name='angletype2')

        ff['A'] = atomtype1
        ff['*'] = atomtype2
        ff['A-B'] = bondtype1
        ff['*-B'] = bondtype2
        ff['A-B-C'] = angletype1
        ff['*-B-C'] = angletype2

        assert ff['Z'] == atomtype2
        assert ff['B-A'] == bondtype1
        assert ff['Z-B'] == bondtype2
        assert ff['C-B-A'] == angletype1
        assert ff['Z-B-C'] == angletype2

    def test_getters_wildcards_conflict(self):
        ff = Forcefield(name='TestFF')
        bondtype1 = BondType(name='bondtype1')
        bondtype2 = BondType(name='bondtype2')

        ff['*-A'] = bondtype1
        ff['*-B'] = bondtype2
        ff['*-C'] = bondtype2
        with pytest.raises(TopologyError):
            ff['A-B']

        assert ff['B-C'] == bondtype2
