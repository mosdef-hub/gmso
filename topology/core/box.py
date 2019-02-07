import numpy as np


def _validate_lengths(lengths):
    lengths = np.asarray(lengths, dtype=float, order='C')
    np.reshape(lengths, newshape=(3,), order='C')

    if np.any(np.less(lengths, [0, 0, 0], )):
        raise ValueError('Negative or 0 value lengths passed.'
                         'Lengths must be a value greater than 0.0'
                         'You passed {}'.format(lengths))
    return lengths


def _validate_angles(angles):
    if angles is None:
        angles = np.asarray([90, 90, 90], dtype=float, order='C')
    else:
        angles = np.asarray(angles, dtype=float, order='C')
        np.reshape(angles, newshape=(3, 1), order='C')
    return angles



class Box(object):
    """A box that bounds a `Topology`.

    Parameters
    ----------
    lengths : array-like, shape(3,), dtype=float
        Lengths of the box [a, b, c].
    angles : array-link, optional, shape(3,), dtype=float
        Interplanar angles, [alpha, beta, gamma], that describe the box shape.


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

    def unit_vectors_from_angles(self):
        (alpha, beta, gamma) = self._angles

        radian_conversion = np.pi / 180.0
        cosa = np.cos(alpha * radian_conversion)
        cosb = np.cos(beta * radian_conversion)
        sinb = np.sin(beta * radian_conversion)
        cosg = np.cos(gamma * radian_conversion)
        sing = np.sin(gamma * radian_conversion)
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

        return np.asarray(box_vec, dtype=np.float)

    def full_vectors_from_angles(self):
        return (self._lengths * self.unit_vectors_from_angles().T).T

    def __repr__(self):
        return "Box(a={}, b={}, c={}, alpha={}, beta={}, gamma={})"\
            .format(*self._lengths, *self._angles)