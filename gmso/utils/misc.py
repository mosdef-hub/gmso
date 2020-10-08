from functools import lru_cache
import unyt as u
from unyt.exceptions import UnitConversionError


#UNYT_CACHE = dict()

def unyt_to_hashable(unyt_or_unyt_iter):
    """Convert a (list of) unyt array or quantity to a hashable tuple."""
    if unyt_or_unyt_iter is None:
        return unyt_or_unyt_iter

    if isinstance(unyt_or_unyt_iter, list):
        hash_coll = tuple(_unyt_to_hashable_single(val) for val in unyt_or_unyt_iter)
        return hash_coll
    else:
        return _unyt_to_hashable_single(unyt_or_unyt_iter)


def _unyt_to_hashable_single(val):
    if isinstance(val, u.unyt_quantity):
        return val.value * conversion_factor(val.units)
    else:
        return None

@lru_cache(maxsize=128)
def conversion_factor(unit):
    return unit.base_value

def ensure_valid_dimensions(quantity_1: u.unyt_quantity,
                            quantity_2: u.unyt_quantity) -> None:
    """Verify unyt quantities are of same dimensions

    This utility will verify that the two unyt quantities provided simplify
    to the same dimensionality. This ensures that each quantity can be
    converted to the desired unit system.

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
