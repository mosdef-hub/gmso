import numpy as np
import pytest
import unyt as u
from unyt.testing import assert_allclose_units

from gmso import Topology
from gmso.core.atom import Atom
from gmso.core.box import Box
from gmso.external.convert_parmed import from_parmed
from gmso.formats.gro import _prepare_atoms
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path
from gmso.utils.io import get_fn, has_mbuild, has_parmed, import_
from gmso.formats.itp import read_itp
from gmso.tests.utils import get_path

class Testitp(BaseTest):
    def test_itp(self):
        #top = Topology.load(get_fn("acn.gro"))
        top=read_itp(get_path("LIQ.itp"))
        assert top != None
        
