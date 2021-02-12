from functools import lru_cache
import unyt as u
from unyt.exceptions import UnitConversionError


def unyt_to_hashable(unyt_or_unyt_iter):
    """Convert a (list of) unyt array or quantity to a hashable tuple."""
    if unyt_or_unyt_iter is None:
        return unyt_or_unyt_iter

    if isinstance(unyt_or_unyt_iter, list):
        hashes = tuple(_unyt_to_hashable_single(val) for val in unyt_or_unyt_iter)
        return hashes
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


def mask_with(iterable, window_size=1, mask='*'):
    """Mask an iterable with the `mask` in a sliding window of size `window_size`

    This method masks an iterable elements with a mask object in a sliding window

    Parameters
    ----------
    iterable: Iterable
        The iterable to mask with
    window_size: int, default=1
        The window size for the mask to be applied
    mask: Any, default='*'
        The mask to apply

    Examples
    --------
        >>> from gmso.utils.misc import mask_with
        >>> list(mask_with(['Ar', 'Ar'], 2))
        [['*', 'Ar'], ['Ar', '*']]
        >>> for masked_list in mask_with(['Ar', 'Xe', 'Xm', 'CH'], 2, mask='_'):
        ...     print('~'.join(masked_list))
        _~_~Xm~CH
        Ar~_~_~CH
        Ar~Xe~_~_

    Yields
    ------
    list
        The masked iterable
    """
    input_list = list(iterable)
    idx = 0

    while idx + window_size <= len(input_list):
        yield input_list[0:idx] +\
              [mask] * window_size +\
              input_list[idx + window_size:]

        idx += 1
