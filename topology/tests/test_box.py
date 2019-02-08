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

    def test_dtype(self):
        box = Box(lengths=np.zeros(3))
        assert box.lengths.dtype == float
        assert box.angles.dtype == float

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

    @pytest.mark.parametrize('lengths, angles', [([3, 3, 3],
        [90.0, 90.0, 90.0]), ([4, 4, 4], [30.0, 60.0, 70.0]),
        ([4, 6, 4], [45.0, 45.0, 75.0])])
    def test_setters_with_lists(self, lengths, angles):
        box = Box(lengths=np.ones(3))
        box.lengths = lengths
        box.angles = angles
        assert (box.lengths == lengths).all()
        assert (box.angles == angles).all()

    def test_default_angles(self):
        box = Box(lengths=np.zeros(3))
        assert (box.angles == np.array([90.0, 90.0, 90.0])).all()

    def test_vectors(self):
        box = Box(lengths=np.ones(3), angles=[40.0, 50.0, 60.0])
        vectors = box.full_vectors_from_angles()
        test_vectors = np.array([[1, 0, 0],
                                [0.5, 0.86603, 0],
                                [0.64278, 0.51344, 0.56852]])
        assert np.isclose(box.vectors(), test_vectors)

    def test_negative_z_vector(self):
        box = Box(lengths=np.ones(3), angles=[-90.0, -40.0, 20.0])
        with pytest.raises(Warning):
            box.unit_vectors_from_angles()