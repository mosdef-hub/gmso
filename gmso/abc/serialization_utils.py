from warnings import warn
from typing import Union
import numpy as np
import unyt as u

__all__ = [
    'unyt_to_dict',
    'dict_to_unyt',
    'GMSOJSONHandler'
]


def unyt_to_dict(unyt_qt: Union[u.unyt_array, u.unyt_quantity]) -> dict:
    if not isinstance(unyt_qt, u.unyt_array):
        raise ValueError(
            'Please provide a value of type unyt array or unyt quantity'
        )
    else:
        numpy_array = unyt_qt.value
        unit = str(unyt_qt.units)
        return {
            'array': numpy_array.tolist(),
            'unit': unit
        }


def dict_to_unyt(dict_obj):
    for key, value in dict_obj.items():
        if value and 'array' in value and 'unit' in value:
            np_array = np.array(value['array'], dtype=float)
            if np_array.shape == tuple():
                unyt_func = u.unyt_quantity
            else:
                unyt_func = u.unyt_array

            dict_obj[key] = unyt_func(
                np_array,
                value['unit']
            )
        elif isinstance(value, dict):
            dict_to_unyt(value)


class JSONHandler:
    def __init__(self):
        self.json_encoders = {}

    def register(self, type_, callable_, override=False):
        """Register a new JSON encoder for GMSO BaseClasses"""
        if type_ not in self.json_encoders:
            self.json_encoders[type_] = callable_
        else:
            if override:
                warn(
                    f'Overriding json serializer for {type_}'
                )
                self.json_encoders[type_] = callable_


GMSOJSONHandler = JSONHandler()
GMSOJSONHandler.register(u.unyt_array, unyt_to_dict)
