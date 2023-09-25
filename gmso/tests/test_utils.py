import pytest
import unyt as u

from gmso.core.atom import Atom
from gmso.core.dihedral import Dihedral
from gmso.utils.io import run_from_ipython
from gmso.utils.misc import unyt_to_hashable
from gmso.utils.sorting import sort_connection_members, sort_connection_strings


def test_unyt_to_hashable():
    hash(unyt_to_hashable(None))
    hash(unyt_to_hashable(1 * u.nm))
    hash(unyt_to_hashable([4, 4] * u.nm))

    assert hash(unyt_to_hashable(1 * u.nm)) == hash(
        unyt_to_hashable(10 * u.angstrom)
    )
    assert hash(unyt_to_hashable(1 * u.kg)) == hash(
        unyt_to_hashable(1000 * u.g)
    )

    assert hash(unyt_to_hashable(1 * u.nm)) != hash(
        unyt_to_hashable(1.01 * u.nm)
    )
    assert hash(unyt_to_hashable(1 * u.nm)) != hash(
        unyt_to_hashable(1.01 * u.second)
    )
    assert hash(unyt_to_hashable(1 * u.nm)) != hash(
        unyt_to_hashable([1, 1] * u.nm)
    )


def test_has_ipython():
    __IPYTHON__ = None
    assert run_from_ipython() is False


def test_sorting():
    with pytest.raises(TypeError):
        sort_connection_members([1])
    atom1 = Atom(name="atom1")
    atom2 = Atom(name="atom2")
    atom3 = Atom(name="atom3")
    atom4 = Atom(name="atom4")

    connect = Dihedral(connection_members=[atom1, atom2, atom3, atom4])
    with pytest.raises(ValueError):
        sort_connection_members(connect, sort_by="error")

    bondList = [atom3.name, atom2.name]
    angleList = [atom3.name, atom2.name, atom1.name]
    dihList = [atom3.name, atom2.name, atom1.name, atom4.name]

    assert sort_connection_strings(bondList) == ("atom2", "atom3")
    assert sort_connection_strings(angleList) == ("atom1", "atom2", "atom3")
    assert sort_connection_strings(dihList) == (
        "atom4",
        "atom1",
        "atom2",
        "atom3",
    )
    assert sort_connection_strings(dihList, improperBool=True) == (
        "atom3",
        "atom1",
        "atom2",
        "atom4",
    )
