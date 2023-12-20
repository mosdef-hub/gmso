from typing import Union
from warnings import warn

import numpy as np
import unyt as u

from gmso.utils.units import GMSO_UnitRegistry

__all__ = ["unyt_to_dict", "dict_to_unyt", "GMSOJSONHandler"]

uregistry = GMSO_UnitRegistry()


def unyt_to_dict(unyt_qt: Union[u.unyt_array, u.unyt_quantity]) -> dict:
    """Convert a unyt quantity into json serializable dictionary"""
    if not isinstance(unyt_qt, u.unyt_array):
        raise TypeError(
            "Please provide a value of type unyt array or unyt quantity"
        )
    else:
        numpy_array = unyt_qt.value
        unit = str(unyt_qt.units)
        return {"array": numpy_array.tolist(), "unit": unit}


def dict_to_unyt(dict_obj) -> None:
    """Recursively convert values of dictionary containing units information to unyt quantities"""
    for key, value in dict_obj.items():
        if isinstance(value, dict):
            if "array" not in value and "unit" not in value:
                dict_to_unyt(value)
            else:
                np_array = np.array(value["array"], dtype=float)
                if np_array.shape == tuple():
                    unyt_func = u.unyt_quantity
                else:
                    unyt_func = u.unyt_array

                dict_obj[key] = unyt_func(
                    np_array, value["unit"], registry=uregistry.reg
                )
