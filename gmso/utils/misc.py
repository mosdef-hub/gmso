import unyt as u
from unyt.exceptions import UnitConversionError


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
        if val.units == u.elementary_charge and val.value == 0:
            val.convert_to_units(u.coulomb)
        val = u.unyt_array(val.reshape(-1))
        return tuple(val.in_base().value)
    else:
        return None


def ensure_valid_dimensions(quantity_1: u.unyt_quantity,
                            quantity_2: u.unyt_quantity) -> None:
    """Ensure that the quantities provided as arguments have valid dimensions

    Parameters
    ----------
    quantity_1: u.unyt_quantity
    quantity_2: u.unyt_quantity
    """
    if quantity_1.units.dimensions != quantity_2.units.dimensions:
        raise UnitConversionError(
            quantity_1.units,
            quantity_1.units.dimensions,
            quantity_2.units,
            quantity_2.units.dimensions
        )
