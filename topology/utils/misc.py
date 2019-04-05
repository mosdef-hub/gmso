import unyt as u


def unyt_to_hashable(val):
    """Convert a unyt array or quantity to a hashable tuple."""
    if val is None:
        return val
    if val.value == 0:
        return tuple((val.reshape(-1).value))
    if isinstance(val, u.unyt_quantity):
        val = u.unyt_array(val.reshape(-1))
    return tuple(val.in_base().value)
