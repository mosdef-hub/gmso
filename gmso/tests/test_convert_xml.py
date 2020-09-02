import numpy as np
import pytest
import unyt as u

from sympy import sympify
from gmso.external.convert_foyer import from_foyer
from gmso.tests.utils import get_path
from gmso.exceptions import ForceFieldParseError
from gmso.tests.base_test import BaseTest
from gmso.utils.io import has_foyer
from gmso.core.forcefield import ForceField

if has_foyer:
   import foyer
   from foyer.tests.utils import get_fn

class TestXMLConversion(BaseTest):

   @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
   @pytest.mark.parametrize("ff", ["fullerene.xml", "oplsaa-periodic.xml", "lj.xml"])
   def test_from_foyer(self, ff):
       from_foyer(get_fn(ff))

   @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
   def test_foyer_version(self, foyer_fullerene):
       assert foyer_fullerene.version == '0.0.1'

   @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
   def test_foyer_scaling(self, foyer_fullerene):
       assert foyer_fullerene.scaling_factors['nonBonded14Scale'] == 1.0
       assert foyer_fullerene.scaling_factors['electrostatics14Scale'] == 1.0

   @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
   def test_foyer_atomtypes(self, foyer_fullerene):
       assert len(foyer_fullerene.atom_types) == 1
       assert 'C' in foyer_fullerene.atom_types

       assert sympify('r') in foyer_fullerene.atom_types['C'].independent_variables
       assert foyer_fullerene.atom_types['C'].parameters['sigma'] == u.unyt_quantity(0.1, u.nm)
       assert foyer_fullerene.atom_types['C'].parameters['ep'] == u.unyt_quantity(0.1, u.kJ / u.mol)
       assert foyer_fullerene.atom_types['C'].mass == u.unyt_quantity(12.01, u.amu)
       assert foyer_fullerene.atom_types['C'].charge == u.unyt_quantity(0.0, u.coulomb)

   def test_empty_foyer_atomtype(self):
       with pytest.raises(ForceFieldParseError):
           from_foyer(get_path("empty_foyer.xml"))
