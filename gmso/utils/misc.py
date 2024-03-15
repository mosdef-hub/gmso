"""Miscellaneous helper methods for GMSO."""

from functools import lru_cache

import unyt as u
from unyt.exceptions import UnitConversionError


@lru_cache(maxsize=128)
def unyt_compare(units1, units2):
    """Check if two parameter values have the same units."""
    units1 = list(units1)
    units2 = list(units2)
    for unit1, unit2 in zip(units1, units2):
        try:
            u.testing.assert_allclose_units(unit1, unit2, rtol=1e-5, atol=1e-8)
        except AssertionError:
            return False
    return True


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
    if quantity_1.units == u.dimensionless:
        return
    elif quantity_1.units.dimensions != quantity_2.units.dimensions:
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


def get_xml_representation(value):
    """Given a value, get its XML representation."""
    if isinstance(value, u.unyt_quantity):
        return str(value.value)
    elif isinstance(value, set):
        return ",".join(value)
    else:
        return str(value)
