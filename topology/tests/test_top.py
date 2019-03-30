import parmed as pmd

from topology.formats.top import write_top
from topology.external.convert_parmed import from_parmed
from topology.tests.base_test import BaseTest
from topology.utils.io import get_fn


class TestTop(BaseTest):
    def test_write_top(self):
        top = from_parmed(pmd.load_file(get_fn('ethane.top'),
                                        xyz=get_fn('ethane.gro')))

        write_top(top, '/Users/mwt/tmp/out.top')
