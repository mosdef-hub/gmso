import warnings

import numpy as np
import unyt as u
from topology.testing.utils import allclose


def _validate_lengths(lengths):
    if not isinstance(lengths, u.unyt_array):
        warnings.warn('Lengths are assumed to be in nm')
        lengths *= u.nm

    input_unit = lengths.units

    lengths = np.asarray(lengths, dtype=float, order='C')
    np.reshape(lengths, newshape=(3,), order='C')

    lengths *= input_unit
    lengths.convert_to_units(u.nm)

    if np.any(np.less(lengths, [0, 0, 0], )):
        raise ValueError('Negative length(s) passed. Lengths must be a value '
                         'greater than 0.0. You passed {}'.format(lengths))

    if np.any(np.equal(lengths, [0, 0, 0], )):
        if lengths[0] > 0 and lengths[1] > 0:
            warnings.warn('A c value of 0 was passed. This will be '
                          'interpreted as a 2-D box.')
        else:
            raise ValueError('Length(s) of value 0 were passed. Lengths must '
                             'be a value greater than 0.0. You passed '
                             '{}'.format(lengths))
    return lengths


def _validate_angles(angles):
    if angles is None:
        angles = np.asarray([90, 90, 90], dtype=float, order='C')
        angles *= u.degree
    else:
        if not isinstance(angles, u.unyt_array):
            warnings.warn('Angles are assumed to be in degrees')
            angles *= u.degree

        input_unit = angles.units

        angles = np.asarray(angles, dtype=float, order='C')
        np.reshape(angles, newshape=(3, 1), order='C')

        angles *= input_unit
        angles.convert_to_units(u.degree)

    return angles



class Box(object):
    """A box that bounds a `Topology`.

    Parameters
    ----------
    lengths : array-like, shape(3,), dtype=float
        Lengths of the box [a, b, c]. Units are assumed to be in nm; if
        passed in as a `unyt_array` it will be converted to nm; if passed in
        as floats, nm is assumed.
    angles : array-like, optional, shape(3,), dtype=float
        Interplanar angles, [alpha, beta, gamma], that describe the box shape.
        Units are assumed to be in degrees; if passed in as a `unyt_array` it
        will be converted to degrees; if passed in as floats, degrees is assumed.


    Attributes
    ----------
    a, b, c : float
        Lengths of the box.
    alpha, beta, gamma : float
        Interplanar angles fully describing the `Box`.

    Methods
    -------
    vectors()
        Output the unit vectors describing the shape of the `Box`.

    """

    def __init__(self, lengths, angles=None):
        """Constructs a `Box`."""

        self._lengths = _validate_lengths(lengths)
        self._angles = _validate_angles(angles)

    @property
    def lengths(self):
        return self._lengths

    @property
    def angles(self):
        return self._angles

    @lengths.setter
    def lengths(self, lengths):
        self._lengths = _validate_lengths(lengths)

    @angles.setter
    def angles(self, angles):
        self._angles = _validate_angles(angles)

    def _unit_vectors_from_angles(self):
        (alpha, beta, gamma) = self.angles

        cosa = np.cos(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)
        cosg = np.cos(gamma)
        sing = np.sin(gamma)
        mat_coef_y = (cosa - cosb * cosg) / sing
        mat_coef_z = np.power(sinb, 2, dtype=float) - \
            np.power(mat_coef_y, 2, dtype=float)

        if mat_coef_z > 0.:
            mat_coef_z = np.sqrt(mat_coef_z)
        else:
            raise Warning('Non-positive z-vector. Angles {} '
                          'do not generate a box with the z-vector in the'
                          'positive z direction'.format(self._angles))

        box_vec = [[1, 0, 0],
                   [cosg, sing, 0],
                   [cosb, mat_coef_y, mat_coef_z]]

        return u.unyt_array(box_vec, u.dimensionless, dtype=np.float)

    def get_vectors(self):
        """ Return the vectors of the box."""
        return (self._lengths * self.get_unit_vectors().T).T
    
    def get_unit_vectors(self):
        """ Return the normalized vectors of the box."""
        return self._unit_vectors_from_angles()

    def __repr__(self):
        return "Box(a={}, b={}, c={}, alpha={}, beta={}, gamma={})"\
            .format(*self._lengths, *self._angles)

    def __eq__(self, other):
        """Compare two boxes for equivalence."""

        if self is other:
            return True

        if not isinstance(other, Box):
            return False

        for not allclose(box.lengths, other.lengths):
                return False

        for not allclose(box.angles, other.angles):
            return False

        return True
