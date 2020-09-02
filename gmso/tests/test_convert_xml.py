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
        assert foyer_fullerene.version == "0.0.1"

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_foyer_scaling(self, foyer_fullerene):
        assert foyer_fullerene.scaling_factors["nonBonded14Scale"] == 1.0
        assert foyer_fullerene.scaling_factors["electrostatics14Scale"] == 1.0

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_foyer_atomtypes(self, foyer_fullerene):
        assert len(foyer_fullerene.atom_types) == 1
        assert "C" in foyer_fullerene.atom_types

        assert sympify("r") in foyer_fullerene.atom_types["C"].independent_variables
        assert foyer_fullerene.atom_types["C"].parameters["sigma"] == u.unyt_quantity(
            0.1, u.nm
        )
        assert foyer_fullerene.atom_types["C"].parameters["ep"] == u.unyt_quantity(
            0.1, u.kJ / u.mol
        )
        assert foyer_fullerene.atom_types["C"].mass == u.unyt_quantity(12.01, u.amu)
        assert foyer_fullerene.atom_types["C"].charge == u.unyt_quantity(0.0, u.coulomb)
        assert foyer_fullerene.atom_types["C"].description == "carbon"
        assert foyer_fullerene.atom_types["C"].definition == "[C;r5;r6]"
        assert foyer_fullerene.atom_types["C"].expression == sympify(
            "ep*(-sigma**6/r**6 + sigma**12/r**12) + q/(e0*r)"
        )

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_foyer_bonds(self, foyer_fullerene):
        assert len(foyer_fullerene.bond_types) == 1
        assert "C~C" in foyer_fullerene.bond_types

        assert sympify("r") in foyer_fullerene.bond_types["C~C"].independent_variables
        assert foyer_fullerene.bond_types["C~C"].parameters["r_eq"] == u.unyt_quantity(
            0.1, u.nm
        )
        assert foyer_fullerene.bond_types["C~C"].parameters["k"] == u.unyt_quantity(
            1000, u.kJ / u.nm ** 2
        )
        assert foyer_fullerene.bond_types["C~C"].member_types == ["C", "C"]

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_foyer_angles(self, foyer_fullerene):
        assert len(foyer_fullerene.angle_types) == 1
        assert "C~C~C" in foyer_fullerene.angle_types

        assert (
            sympify("theta")
            in foyer_fullerene.angle_types["C~C~C"].independent_variables
        )
        assert foyer_fullerene.angle_types["C~C~C"].parameters["k"] == u.unyt_quantity(
            1000, u.kJ / u.rad ** 2
        )
        assert foyer_fullerene.angle_types["C~C~C"].parameters[
            "theta_eq"
        ] == u.unyt_quantity(3.141592, u.rad)

    def test_empty_foyer_atomtype(self):
        with pytest.raises(ForceFieldParseError):
            from_foyer(get_path("empty_foyer.xml"))
