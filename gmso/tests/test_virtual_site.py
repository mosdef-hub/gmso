import pytest

import mbuild as mb
import numpy as np
import unyt as u
from pydantic import ValidationError

from gmso.core.virtual_site import VirtualSite 
from gmso.core.atom import Atom
from gmso.tests.base_test import BaseTest


class TestVirturalSite(BaseTest):
    def test_new_site(self, water_system):
        v_site = VirtualSite(parent_atoms = water_system.sites)
        assert len(v_site.parent_atoms) == 3
        for site in v_site.parent_atoms:
            assert site in water_system.sites
