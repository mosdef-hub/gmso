import pytest
import numpy as np
from topology.core.box import Box


class TestBox():

    def test_init_lengths(self):
        box = Box(lengths=np.ones(3))
        assert np.array_equal(box.lengths, np.ones(3))

    def test_init_angles(self):
        box = Box(lengths=np.ones(3), angles=[40.0, 50.0, 60.0])
        assert np.array_equal(box.angles, [40.0, 50.0, 60.0])

    @pytest.mark.parametrize('lengths', [[4.0, 4.0, 0.0], [4.0, 5.0, -1.0]])
    def test_bad_lengths(self, lengths):
        with pytest.raises(ValueError):
            box = Box(lengths=lengths, angles=90*np.ones(3))

    def test_dtype(self):
        box = Box(lengths=np.ones(3))
        assert box.lengths.dtype == float

    def test_lengths_setter(self):
        box = Box(lengths=np.ones(3))
        box.lengths = 2 * np.ones(3)
        assert (box.lengths == 2 * np.ones(3)).all()

    @pytest.mark.parametrize('angles', [[40.0, 50.0, 60.0],
        [30.0, 60.0, 70.0], [45.0, 45.0, 75.0]])
    def test_angles_setter(self, angles):
        box = Box(lengths=np.ones(3), angles=90*np.ones(3))
        box.angles = angles
        assert (box.angles == angles).all()

    @pytest.mark.parametrize('lengths', [[3, 3, 3], [4, 4, 4],
        [4, 6, 4]])
    def test_setters_with_lists(self, lengths):
        box = Box(lengths=np.ones(3))
        box.lengths = lengths
        assert (box.lengths == lengths).all()

    def test_default_angles(self):
        box = Box(lengths=np.ones(3))
        assert (box.angles == np.array([90.0, 90.0, 90.0])).all()

    def test_vectors(self):
        box = Box(lengths=np.ones(3), angles=[40.0, 50.0, 60.0])
        vectors = box.full_vectors_from_angles()
        test_vectors = np.array([[1, 0, 0],
                                [0.5, 0.86603, 0],
                                [0.64278, 0.51344, 0.56852]])
        assert np.isclose(box.vectors(), test_vectors)
