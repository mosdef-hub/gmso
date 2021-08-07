from copy import deepcopy

import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

from gmso.core.box import Box
from gmso.tests.base_test import BaseTest


class TestBox(BaseTest):
    def test_init_lengths(self, lengths):
        box = Box(lengths=lengths)
        assert np.array_equal(box.lengths, lengths)

    def test_init_angles(self, lengths, angles):
        box = Box(lengths=lengths, angles=angles)
        assert np.array_equal(box.angles, angles)

    @pytest.mark.parametrize("lengths", [[0.0, 4.0, 4.0], [4.0, 5.0, -1.0]])
    def test_bad_lengths(self, lengths, angles):
        lengths *= u.nm
        with pytest.raises(ValueError):
            Box(lengths=lengths, angles=angles)

    def test_build_2D_Box(self):
        with pytest.warns(UserWarning):
            Box(lengths=u.nm * [4, 4, 0])

    def test_dtype(self, box):
        assert box.lengths.dtype == float
        assert isinstance(box.lengths, u.unyt_array)
        assert isinstance(box.lengths, np.ndarray)
        assert box.angles.dtype == float
        assert isinstance(box.angles, u.unyt_array)
        assert isinstance(box.angles, np.ndarray)

    def test_lengths_setter(self, box, lengths):
        box.lengths = 2 * u.nm * np.ones(3)
        assert (box.lengths == 2 * u.nm * np.ones(3)).all()

    @pytest.mark.parametrize(
        "angles", [[40.0, 50.0, 60.0], [30.0, 60.0, 70.0], [45.0, 45.0, 75.0]]
    )
    def test_angles_setter(self, lengths, angles):
        box = Box(lengths=lengths, angles=u.degree * np.ones(3))
        angles *= u.degree
        box.angles = angles
        assert (box.angles == angles).all()

    @pytest.mark.parametrize("lengths", [[3, 3, 3], [4, 4, 4], [4, 6, 4]])
    def test_setters_with_lists(self, lengths):
        box = Box(lengths=u.nm * np.ones(3))
        lengths *= u.nm
        box.lengths = lengths
        assert (box.lengths == lengths).all()

    def test_default_angles(self, box):
        assert (box.angles == u.degree * np.array([90.0, 90.0, 90.0])).all()

    def test_unit_vectors(self):
        box = Box(
            lengths=u.nm * np.ones(3), angles=u.degree * [40.0, 50.0, 60.0]
        )
        vectors = box.get_unit_vectors()
        assert vectors.units.is_dimensionless

    def test_scaling_unit_vectors(self):
        box = Box(
            lengths=u.unyt_array((2, 2, 2), u.nm),
            angles=u.degree * [40.0, 50.0, 60.0],
        )
        u_vectors = box.get_unit_vectors()
        scaled_u_vectors = (u_vectors.T * box.lengths).T

        scaled_vectors = box.get_vectors()

        assert_allclose_units(
            scaled_vectors, scaled_u_vectors, rtol=1e-5, atol=u.nm * 1e-3
        )
        assert scaled_u_vectors.units == u.nm

    def test_scaled_vectors(self):
        box = Box(
            lengths=u.unyt_array((2, 2, 2), u.nm),
            angles=u.degree * [40.0, 50.0, 60.0],
        )
        vectors = box.get_vectors()
        test_vectors = np.array(
            [[1, 0, 0], [0.5, 0.86603, 0], [0.64278, 0.51344, 0.56852]]
        )
        test_vectors = (test_vectors.T * box.lengths).T
        assert_allclose_units(
            vectors, test_vectors, rtol=1e-5, atol=u.nm * 1e-3
        )
        assert vectors.units == u.nm

    def test_eq(self, box):
        assert box == box

    def test_eq_bad_lengths(self, box):
        diff_lengths = deepcopy(box)
        diff_lengths.lengths = u.nm * [5.0, 5.0, 5.0]
        assert box != diff_lengths

    def test_eq_bad_angles(self, box):
        diff_angles = deepcopy(box)
        diff_angles.angles = u.degree * [90, 90, 120]
        assert box != diff_angles

    def test_list_to_unyt_array(self):
        length = 5 * u.nm
        box = Box([length, length, length])

        assert isinstance(box.lengths, u.unyt_array)
