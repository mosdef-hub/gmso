from gmso.formats.itp import read_itp
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path


class Testitp(BaseTest):
    def test_itp(self):
        # top = Topology.load(get_fn("acn.gro"))
        top = read_itp(get_path("LIQ.itp"))
        assert top != None
