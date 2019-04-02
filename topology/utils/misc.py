import unyt as u


def unyt_to_hashable(val):
    """Convert a unyt array or quantity to a hashable tuple."""
    if isinstance(val, u.unyt_quantity):
        val = u.unyt_array(val.reshape(-1))
    return hash(tuple(val.in_base().value))
