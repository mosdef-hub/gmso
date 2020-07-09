import numpy as np
import pytest

from gmso.external.convert_foyer import from_foyer, to_foyer
from gmso.tests.utils import get_path
from gmso.exceptions import ForceFieldParseError
from gmso.tests.base_test import BaseTest
from gmso.utils.io import has_foyer

if has_foyer:
    import foyer
    from foyer.tests.utils import get_fn

class TestXMLConversion(BaseTest):

    @pytest.mark.parametrize("ff", ["spce.xml", "carbon.xml", "ethylene.xml"])
    def test_to_foyer(self, ff):
        to_foyer(get_path(ff))

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    @pytest.mark.parametrize("ff", ["fullerene.xml", "oplsaa-periodic.xml", "lj.xml"])
    def test_from_foyer(self, ff):
        from_foyer(get_fn(ff))

    def test_empty_foyer_atomtype(self):
        with pytest.raises(ForceFieldParseError):
            from_foyer(get_path("empty.xml"))
