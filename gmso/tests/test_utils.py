import unyt as u

from gmso.utils.io import run_from_ipython
from gmso.utils.misc import are_equal_unyt_dicts, unyt_to_hashable


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


def test_are_equal_unyt_dicts():
    u1 = {"a": 2.0 * u.nm, "b": 3.5 * u.nm}
    u2 = {"c": 2.0 * u.nm, "d": 3.5 * u.nm}
    assert are_equal_unyt_dicts(u1, u2) is False
