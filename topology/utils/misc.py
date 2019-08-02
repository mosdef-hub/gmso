import unyt as u


def unyt_to_hashable(val):
    """Convert a unyt array or quantity to a hashable tuple."""
    if val is None:
        return val
    # TODO: More elegantly handle hashing charges, see #123
    if isinstance(val, u.unyt_quantity):
        if val.units == u.elementary_charge and u.value == 0:
            val.convert_to_units(u.coulomb)
        val = u.unyt_array(val.reshape(-1))
    return tuple(val.in_base().value)
