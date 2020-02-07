import unyt as u


def unyt_to_hashable(unyt_or_unyt_iter):
    """Convert a (list of) unyt array or quantity to a hashable tuple."""
    if unyt_or_unyt_iter is None:
        return unyt_or_unyt_iter
    hash_coll = []
    if isinstance(unyt_or_unyt_iter, list):
        for val in unyt_or_unyt_iter:
            hash_coll.append(_unyt_to_hashable_single(val))
        return tuple(hash_coll)
    else:
        return _unyt_to_hashable_single(unyt_or_unyt_iter)


def _unyt_to_hashable_single(val):
    if isinstance(val, u.unyt_quantity):
        if val.units == u.elementary_charge and u.value == 0:
            val.convert_to_units(u.coulomb)
        val = u.unyt_array(val.reshape(-1))
        return tuple(val.in_base().value)
    else:
        return None
