import unyt as u

from topology.utils.misc import unyt_to_hashable

def test_unyt_to_hashable():
    hash(unyt_to_hashable(None))
    hash(unyt_to_hashable(1 * u.nm))
    hash(unyt_to_hashable([4, 4] * u.nm))

    assert hash(unyt_to_hashable(1 * u.nm)) == hash(unyt_to_hashable(10 * u.angstrom))
    assert hash(unyt_to_hashable(1 * u.kg)) == hash(unyt_to_hashable(1000 * u.g))

    assert hash(unyt_to_hashable(1 * u.nm)) != hash(unyt_to_hashable(1.01 * u.nm))
    assert hash(unyt_to_hashable(1 * u.nm)) != hash(unyt_to_hashable(1.01 * u.second))
    assert hash(unyt_to_hashable(1 * u.nm)) != hash(unyt_to_hashable([1, 1] * u.nm))
