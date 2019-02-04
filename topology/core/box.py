import numpy as np


class Box(object):
    """A box representing the bounds of the topology.

    Attributes
    ----------
    mins : np.ndarray, shape=(3,), dtype=float
        Minimum x, y, z coordinates.
    maxs : np.ndarray, shape=(3,), dtype=float
        Maximum x, y, z coordinates.
    lengths : np.ndarray, shape(3,), dtype=float
        Box length in x, y and z directions.
    vectors : np.ndarray, shape(3,3), dtype=float
        Unit vectors that define the edges of the box.
    positions : np.ndarray, shape(n_sites, 3), dtype=float
        XYZ coordinates of the sites in the topology.

    """
    def __init__(self, lengths=None, mins=None, maxs=None,
                 angles=None, vectors=None, positions=None):
        lengths, mins, maxs, angles, vectors = self._validate(
            lengths=lengths, mins=mins, maxs=maxs,
            angles=angles, vectors=vectors, positions=positions)
        self._lengths = lengths
        self._mins = mins
        self._maxs = maxs
        self._angles = angles
        self._vectors = vectors
        self._positions = positions
        if lengths is not None:
            assert mins is None and maxs is None
            self._mins = np.array([0.0, 0.0, 0.0])
            self._maxs = np.array(lengths, dtype=np.float)
            self._lengths = np.array(lengths, dtype=np.float)
        elif maxs is not None:
            assert mins is not None and lengths is None
            self._mins = np.array(mins, dtype=np.float)
            self._maxs = np.array(maxs, dtype=np.float)
            self._lengths = self.maxs - self.mins
        else:
            raise ValueError("Either provide `lengths` or `mins` and `maxs`."
                             "You provided: lengths={} mins={} maxs={}".format(lengths, mins, maxs))
        if angles is None:
            angles = np.array([90.0, 90.0, 90.0])
        elif isinstance(angles, (list, np.ndarray)):
            angles = np.array(angles, dtype=np.float)
        self._angles = angles

    def _validate(self, lengths=None, mins=None, maxs=None,
                  angles=None, vectors=None, positions=None):
        """
        Validate the inputs for the box class.


        :return:
        mins : np.ndarray, shape=(3,), dtype=float
            Minimum x, y, z coordinates.
        maxs : np.ndarray, shape=(3,), dtype=float
            Maximum x, y, z coordinates.
        lengths : np.ndarray, shape(3,), dtype=float
            Box length in x, y and z directions.
        vectors : np.ndarray, shape(3,3), dtype=float
            Unit vectors that define the edges of the box.
        positions : np.ndarray, shape (n_sites, 3), dtype=float

        """

        if lengths is not None:



    @property
    def mins(self):
        return self._mins

    @property
    def maxs(self):
        return self._maxs

    @property
    def lengths(self):
        return self._lengths

    @property
    def angles(self):
        return self._angles

    @property
    def vectors(self):
        return self._vectors

    @property
    def positions(self):
        return self._positions

    @mins.setter
    def mins(self, mins):
        if isinstance(mins, list):
            mins = np.array(mins, dtype=np.float)
        assert mins.shape == (3, )
        self._mins = mins
        self._lengths = self.maxs - self.mins

    @maxs.setter
    def maxs(self, maxes):
        if isinstance(maxes, list):
            maxes = np.array(maxes, dtype=np.float)
        assert maxes.shape == (3, )
        self._maxs = maxes
        self._lengths = self.maxs - self.mins

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
        return "Box(mins={}, maxs={}, angles={})".format(self.mins, self.maxs, self.angles)