import pytest
import numpy as np
import unyt as u

from topology.core.box import Box
from topology.tests.base_test import BaseTest
from topology.testing.utils import allclose


class TestBox(BaseTest):
    def test_init_lengths(self, lengths):
        box = Box(lengths=lengths)
        assert np.array_equal(box.lengths, lengths)

    def test_init_angles(self, lengths, angles):
        box = Box(lengths=lengths, angles=angles)
        assert np.array_equal(box.angles, angles)

    @pytest.mark.parametrize('lengths', [[0.0, 4.0, 4.0], [4.0, 5.0, -1.0]])
    def test_bad_lengths(self, lengths, angles):
        lengths *= u.nm
        with pytest.raises(ValueError):
            box = Box(lengths=lengths, angles=angles)

    def test_build_2D_Box(self):
        with pytest.warns(UserWarning):
            box = Box(lengths=u.nm * [4, 4, 0])

    def test_dtype(self, box):
        assert box.lengths.dtype == float
        assert isinstance(box.lengths, u.unyt_array)
        assert isinstance(box.lengths, np.ndarray)
        assert box.angles.dtype == float
        assert isinstance(box.angles, u.unyt_array)
        assert isinstance(box.angles, np.ndarray)

    def test_lengths_setter(self, box, lengths):
        box.lengths = 2*u.nm*np.ones(3)
        assert (box.lengths == 2*u.nm*np.ones(3)).all()

    @pytest.mark.parametrize('angles', [[40.0, 50.0, 60.0],
        [30.0, 60.0, 70.0], [45.0, 45.0, 75.0]])
    def test_angles_setter(self, lengths, angles):
        box = Box(lengths=lengths, angles=u.degree*np.ones(3))
        angles *= u.degree
        box.angles = angles
        assert (box.angles == angles).all()

    @pytest.mark.parametrize('lengths', [[3, 3, 3], [4, 4, 4],
        [4, 6, 4]])
    def test_setters_with_lists(self, lengths):
        box = Box(lengths=u.nm*np.ones(3))
        lengths *= u.nm
        box.lengths = lengths
        assert (box.lengths == lengths).all()

    def test_default_angles(self, box):
        assert (box.angles == u.degree * np.array([90.0, 90.0, 90.0])).all()

    def test_vectors(self):
        box = Box(lengths=u.nm*np.ones(3), angles=u.degree*[40.0, 50.0, 60.0])
        vectors = box.full_vectors_from_angles()
        test_vectors = np.array([[1, 0, 0],
                                [0.5, 0.86603, 0],
                                [0.64278, 0.51344, 0.56852]])
        test_vectors *= u.nm
        assert allclose(vectors, test_vectors, atol=u.nm*1e-3)
