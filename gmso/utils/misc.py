"""Miscellaneous helper methods for GMSO."""
from functools import lru_cache

import unyt as u
from unyt.exceptions import UnitConversionError


def unyt_to_hashable(unyt_or_unyt_iter):
    """Convert a (list of) unyt array or quantity to a hashable tuple."""
    if unyt_or_unyt_iter is None:
        return unyt_or_unyt_iter

    if isinstance(unyt_or_unyt_iter, list):
        hashes = tuple(
            _unyt_to_hashable_single(val) for val in unyt_or_unyt_iter
        )
        return hashes
    else:
        return _unyt_to_hashable_single(unyt_or_unyt_iter)


def _unyt_to_hashable_single(val):
    """Convert a unit quantity to a hashable value."""
    if isinstance(val, u.unyt_quantity):
        return val.value * conversion_factor(val.units)
    else:
        return None


@lru_cache(maxsize=128)
def conversion_factor(unit):
    """Store the conversion factor for various unit conversions."""
    return unit.base_value


def ensure_valid_dimensions(
    quantity_1: u.unyt_quantity, quantity_2: u.unyt_quantity
) -> None:
    """Verify unyt quantities are of same dimensions.

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
            quantity_2.units.dimensions,
        )


def validate_type(iterator, type_):
    """Validate all the elements of the iterable are of a particular type."""
    for item in iterator:
        if not isinstance(item, type_):
            raise TypeError(
                f"Expected {item} to be of type {type_.__name__} but got"
                f" {type(item).__name__} instead."
            )


def mask_with(iterable, window_size=1, mask="*"):
    """Mask an iterable with the `mask` in a circular sliding window of size `window_size`.

    This method masks an iterable elements with a mask object in a circular sliding window

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
        >>> list(mask_with(['Ar', 'Ar'], 1))
        [['*', 'Ar'], ['Ar', '*']]
        >>> for masked_list in mask_with(['Ar', 'Xe', 'Xm', 'CH'], 2, mask='_'):
        ...     print('~'.join(masked_list))
        _~_~Xm~CH
        Ar~_~_~CH
        Ar~Xe~_~_
        _~Xe~Xm~_

    Yields
    ------
    list
        The masked iterable
    """
    input_list = list(iterable)
    idx = 0
    first = None
    while idx < len(input_list):
        mask_idxes = set(
            (idx + j) % len(input_list) for j in range(window_size)
        )
        to_yield = [
            mask if j in mask_idxes else input_list[j]
            for j in range(len(input_list))
        ]
        if to_yield == first:
            break
        if idx == 0:
            first = to_yield

        idx += 1
        yield to_yield


def are_equal_unyt_dicts(u1, u2):
    """Compare two dictionaries of unyt quantities/arrays.

    This method compares two dictionaries (`u1` and `u2`) of
    `unyt_quantities` and returns True if:
        * u1 and u2 have the exact same key set
        * for each key, the value in u1 and u2 have the same unyt quantity

    Notes
    -----
    Type checks are not performed for the sake of removing unnecessary
    branching and it is incumbent upon the callee to incorporate correct
    messages.
    """
    if u1.keys() != u2.keys():
        return False
    else:
        for k, v in u1.items():
            if not u.allclose_units(v, u2[k]):
                return False

        return True
