import numpy as np


class Box(object):
    """A box representing the bounds of the topology.

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

        lengths = np.asarray(lengths, dtype=float, order='C')
        np.reshape(lengths, newshape=(3,), order='C')

        if np.any(np.less(lengths, [0, 0, 0],)):
            raise ValueError('Negative or 0 value lengths passed.'
                             'Lengths must be a value greater than 0.0'
                             'You passed {}'.format(lengths))

        if angles is None:
            angles = np.asarray([90, 90, 90], dtype=float, order='C')
        else:
            angles = np.asarray(angles, dtype=float, order='C')
            np.reshape(angles, newshape=(3,1), order='C')

        self._lengths = lengths
        self._angles = angles

    @property
    def lengths(self):
        return self._lengths

    @property
    def angles(self):
        return self._angles


    @lengths.setter
    def lengths(self, lengths):
        if isinstance(lengths, list):
            lengths = np.array(lengths, dtype=np.float)
        assert lengths.shape == (3, )
        self._maxs += 0.5*lengths - 0.5*self.lengths
        self._mins -= 0.5*lengths - 0.5*self.lengths
        self._lengths = lengths

    @angles.setter
    def angles(self, angles):
        if isinstance(angles, list):
            angles = np.array(angles, dtype=np.float)
        assert angles.shape == (3, )
        self._angles = angles

    def __repr__(self):
        return "Box(a={}, b={}, c={}, alpha={}, beta={}, gamma={})"\
            .format(*self._lengths, *self._angles)