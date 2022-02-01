from pathlib import Path

import pytest
import unyt as u
from sympy import sympify

from gmso.exceptions import ForceFieldParseError
from gmso.external.convert_foyer_xml import from_foyer_xml
from gmso.tests.base_test import BaseTest
from gmso.tests.utils import get_path
from gmso.utils.io import has_foyer

if has_foyer:
    from foyer.tests.utils import get_fn

parameterized_ffs = ["fullerene.xml", "oplsaa-periodic.xml", "lj.xml"]


@pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
class TestXMLConversion(BaseTest):
    @pytest.mark.parametrize("ff", parameterized_ffs)
    def test_from_foyer(self, ff):
        from_foyer_xml(get_fn(ff), overwrite=True)

    @pytest.mark.parametrize("ff", parameterized_ffs)
    def test_from_foyer_overwrite_false(self, ff):
        from_foyer_xml(get_fn(ff), overwrite=False)
        with pytest.raises(FileExistsError):
            from_foyer_xml(get_fn(ff), overwrite=False)

    @pytest.mark.parametrize("ff", parameterized_ffs)
    def test_from_foyer_different_name(self, ff):
        from_foyer_xml(get_fn(ff), f"{ff}-gmso-converted.xml", overwrite=True)

    @pytest.mark.parametrize("ff", parameterized_ffs)
    def test_from_foyer_validate_foyer(self, ff):
        from_foyer_xml(
            get_fn(ff),
            f"{ff}-gmso-converted.xml",
            overwrite=True,
            validate_foyer=True,
        )

    @pytest.mark.parametrize("ff", parameterized_ffs)
    def test_foyer_pathlib(self, ff):
        file_path = Path(get_fn(ff))
        from_foyer_xml(file_path, overwrite=True)

    def test_foyer_file_not_found(self):
        file_path = "dummy_name.xml"
        with pytest.raises(FileNotFoundError):
            from_foyer_xml(file_path, overwrite=True)

    def test_foyer_version(self, foyer_fullerene):
        assert foyer_fullerene.version == "0.0.1"

    def test_foyer_14scale(self, foyer_fullerene):
        assert foyer_fullerene.scaling_factors["electrostatics14Scale"] == 1.0
        assert foyer_fullerene.scaling_factors["nonBonded14Scale"] == 1.0

    def test_foyer_scaling(self, foyer_fullerene):
        assert foyer_fullerene.scaling_factors["nonBonded14Scale"] == 1.0
        assert foyer_fullerene.scaling_factors["electrostatics14Scale"] == 1.0

    def test_foyer_atomtypes(self, foyer_fullerene):
        assert len(foyer_fullerene.atom_types) == 1
        assert "C" in foyer_fullerene.atom_types

        assert (
            sympify("r")
            in foyer_fullerene.atom_types["C"].independent_variables
        )
        assert foyer_fullerene.atom_types["C"].parameters[
            "sigma"
        ] == u.unyt_quantity(0.1, u.nm)
        assert foyer_fullerene.atom_types["C"].parameters[
            "epsilon"
        ] == u.unyt_quantity(0.1, u.kJ / u.mol)
        assert foyer_fullerene.atom_types["C"].mass == u.unyt_quantity(
            12.01, u.amu
        )
        assert foyer_fullerene.atom_types["C"].charge == u.unyt_quantity(
            0.0, u.coulomb
        )
        assert foyer_fullerene.atom_types["C"].description == "carbon"
        assert foyer_fullerene.atom_types["C"].definition == "[C;r5;r6]"
        assert foyer_fullerene.atom_types["C"].expression == sympify(
            "epsilon*(-sigma**6/r**6 + sigma**12/r**12)"
        )

    def test_foyer_bonds(self, foyer_fullerene):
        assert len(foyer_fullerene.bond_types) == 1
        assert "C~C" in foyer_fullerene.bond_types

        assert (
            sympify("r")
            in foyer_fullerene.bond_types["C~C"].independent_variables
        )
        assert foyer_fullerene.bond_types["C~C"].parameters[
            "r_eq"
        ] == u.unyt_quantity(0.1, u.nm)
        assert foyer_fullerene.bond_types["C~C"].parameters[
            "k"
        ] == u.unyt_quantity(1000, u.kJ / u.mol / u.nm**2)
        assert foyer_fullerene.bond_types["C~C"].member_classes == ("C", "C")

    def test_foyer_angles(self, foyer_fullerene):
        assert len(foyer_fullerene.angle_types) == 1
        assert "C~C~C" in foyer_fullerene.angle_types

        assert (
            sympify("theta")
            in foyer_fullerene.angle_types["C~C~C"].independent_variables
        )
        assert foyer_fullerene.angle_types["C~C~C"].parameters[
            "k"
        ] == u.unyt_quantity(1000, u.kJ / u.mol / u.rad**2)
        assert foyer_fullerene.angle_types["C~C~C"].parameters[
            "theta_eq"
        ] == u.unyt_quantity(3.141592, u.rad)
        assert foyer_fullerene.angle_types["C~C~C"].member_classes == (
            "C",
            "C",
            "C",
        )

    def test_foyer_dihedrals(self, foyer_periodic):
        assert len(foyer_periodic.dihedral_types) == 4
        assert (
            "opls_140~opls_135~opls_135~opls_140"
            in foyer_periodic.dihedral_types
        )

        assert (
            sympify("phi")
            in foyer_periodic.dihedral_types[
                "opls_140~opls_135~opls_135~opls_140"
            ].independent_variables
        )
        assert foyer_periodic.dihedral_types[
            "opls_140~opls_135~opls_135~opls_140"
        ].parameters["k"] == u.unyt_quantity(3.1, u.kJ / u.mol)
        assert foyer_periodic.dihedral_types[
            "opls_140~opls_135~opls_135~opls_140"
        ].parameters["n"] == u.unyt_quantity(1, u.dimensionless)
        assert foyer_periodic.dihedral_types[
            "opls_140~opls_135~opls_135~opls_140"
        ].parameters["delta"] == u.unyt_quantity(3.14, u.rad)
        assert foyer_periodic.dihedral_types[
            "opls_140~opls_135~opls_135~opls_140"
        ].member_types == ("opls_140", "opls_135", "opls_135", "opls_140")

    def test_foyer_urey_bradley(self, foyer_urey_bradley):
        assert foyer_urey_bradley.angle_types["OBL~CL~CTL2"] is not None

    def test_foyer_rb_torsion(self, foyer_rb_torsion):
        assert foyer_rb_torsion.dihedral_types["HC~CT~CT~HC"] is not None

    def test_empty_foyer_atomtype(self):
        with pytest.raises(ForceFieldParseError):
            from_foyer_xml(get_path("empty_foyer.xml"))
