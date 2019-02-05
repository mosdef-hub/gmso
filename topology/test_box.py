import pytest
import numpy as np
from core.box import Box


class TestBox():

    def test_init_lengths(self):
        box = Box(lengths=np.ones(3))
        assert np.array_equal(box.lengths, np.ones(3))
        assert np.array_equal(box.mins, np.zeros(3))
        assert np.array_equal(box.maxs, np.ones(3))

    def test_init_angles(self):
        box = Box(lengths=np.ones(3), angles=[40.0, 50.0, 60.0])
        assert np.array_equal(box.angles, [40.0, 50.0, 60.0])

    def test_dtype(self):
        box = Box(lengths=np.zeros(3)))
        assert box.lengths.dtype == np.float64
        assert box.mins.dtype == np.float64
        assert box.maxs.dtype == np.float64

    def test_mins_setter(self):
        box = Box(mins=np.zeros(3), maxs=2 * np.ones(3))
        box.mins = np.ones(3)
        assert (box.mins == np.ones(3)).all()
        assert (box.maxs - box.mins == np.ones(3)).all()
        assert (box.lengths == np.ones(3)).all()

    def test_maxs_setter(self):
        box = Box(mins=np.zeros(3), maxs=np.ones(3))
        box.maxs = 2 * np.ones(3)
        assert (box.maxs == 2 * np.ones(3)).all()
        assert (box.maxs - box.mins == 2 * np.ones(3)).all()
        assert (box.lengths == 2 * np.ones(3)).all()

    def test_lengths_setter(self):
        box = Box(lengths=np.ones(3))
        box.lengths = 2 * np.ones(3)
        assert (box.lengths == 2 * np.ones(3)).all()
        assert (box.maxs - box.mins == 2 * np.ones(3)).all()

    def test_angles_setter(self):
        box = Box(lengths=np.ones(3), angles=90*np.ones(3))
        box.angles = np.array([60.0, 120.0, 60.0])
        assert (box.angles == np.array([60.0, 120.0, 60.0])).all()

    def test_setters_with_lists(self):
        box = Box(lengths=np.ones(3))
        box.lengths = [4, 4, 4]
        assert (box.lengths == 4 * np.ones(3)).all()

    def test_default_angles(self):
        box = Box(lengths=np.zeros(3))
        assert (box.angles == np.array([90.0, 90.0, 90.0])).all()
