import numpy as np

from topology.utils.decorators import confirm_dict_existence
from topology.exceptions import TopologyError


class TabulatedPotential(object):
    """A class representing an potential represented by tabulated values."""

    def __init__(self, name="Potential", r=None, v=None, topology=None):
        self._name = name
        self._r = r
        self._v = v

        if topology is not None:
            self._topology = topology
        else:
            self._topology = None

    @property
    def name(self):
        return self._name

    @name.setter
    @confirm_dict_existence
    def name(self, val):
        self._name = val

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, val):
        val = _preprocess_array(val)
        _check_consistency(val, self._v)
        self._r = val

    @property
    def v(self):
        return self._v

    @v.setter
    def v(self, val):
        val = _preprocess_array(val)
        _check_consistency(self._r, val)
        self._v = val

    def set_r_and_v(self, r, v):
        r_ = _preprocess_array(r)
        v_ = _preprocess_array(v)
        _check_consistency(r_, v_)
        self._r = r_
        self._v = v_


def _preprocess_array(arr):
    arr = np.asarray(arr)

    if arr.ndim != 1:
        raise TopologyError(
            'Array found to have non-1 number of dimensions after conversion'
            'to NumPy. Number of dimensions found: {arr.ndim}'
        )

    return arr


def _check_consistency(a, b):
    if a.shape != b.shape:
        raise TopologyError(
            'Trying to set an array of length unequal to other attribute. See '
            'TabulatedPotential.set_r_and_v() to set both r and v in one step'
        )
