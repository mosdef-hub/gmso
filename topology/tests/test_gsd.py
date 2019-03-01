import numpy as np
import unyt as u
import parmed as pmd
import pytest

from topology.formats.gsd import write_gsd
from topology.external.convert_parmed import from_parmed
from topology.tests.base_test import BaseTest
from topology.utils.io import get_fn
from topology.testing.utils import allclose


class TestGsd(BaseTest):
    def test_write_gsd(self):
        top = from_parmed(pmd.load_file(get_fn('ethane.top'), 
            xyz=get_fn('ethane.gro')))

        write_gsd(top, 'out.gsd')

    def test_write_gsd_non_orthogonal(self):
        top = from_parmed(pmd.load_file(get_fn('ethane.top'), 
            xyz=get_fn('ethane.gro')))
        top.box.angles = u.degree * [90, 90, 120]

        write_gsd(top, 'out.gsd')
