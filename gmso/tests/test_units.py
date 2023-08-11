import re

import pytest
import unyt as u

from gmso.tests.base_test import BaseTest
from gmso.utils.misc import unyt_to_hashable
from gmso.utils.units import LAMMPS_UnitSystems


class TestUnitHandling(BaseTest):
    @pytest.fixture
    def real_usys(self):
        return LAMMPS_UnitSystems("real")

    def test_unyt_to_hashable(self):
        hash(unyt_to_hashable(None))
        hash(unyt_to_hashable(1 * u.nm))
        hash(unyt_to_hashable([4, 4] * u.nm))

        assert hash(unyt_to_hashable(1 * u.nm)) == hash(
            unyt_to_hashable(10 * u.angstrom)
        )
        assert hash(unyt_to_hashable(1 * u.kg)) == hash(
            unyt_to_hashable(1000 * u.g)
        )

        assert hash(unyt_to_hashable(1 * u.nm)) != hash(
            unyt_to_hashable(1.01 * u.nm)
        )
        assert hash(unyt_to_hashable(1 * u.nm)) != hash(
            unyt_to_hashable(1.01 * u.second)
        )
        assert hash(unyt_to_hashable(1 * u.nm)) != hash(
            unyt_to_hashable([1, 1] * u.nm)
        )

    """
    Utilities to make
    [a] register units needed for unit systems
    [a] need to be able to generate unit systems
    [a] take a unyt and check for energy
    [a] take a unyt and check for electron volts
    [?] take a unyt and look for thermal energy
    [a] get all dimensions from a unit
    [a] convert a unit using a base system
    [a] return a rounded float for a unit
    [] be able to write out a unit without the conversion
    [] attach units to a float with a unit system and dimensions
    [] should have associated somewhere to handle units outside of unit system
    # units should be done as a top level conversion
    # need to implement this module into the writers
    # can probably subclass some functions out of LAMMPS_UnitSystems

    Tests to make

    """

    def test_unit_conversion(self, real_usys):
        parameter = 0.001 * u.Unit("mm")
        n_decimals = 5
        outStr = real_usys.convert_parameter(parameter, n_decimals=n_decimals)
        assert float(outStr) == 10000.00000
        assert outStr[::-1].find(".") == n_decimals

        parameter = 1 * u.Unit("nm")
        n_decimals = 5
        outStr = real_usys.convert_parameter(parameter, n_decimals=n_decimals)
        assert float(outStr) == 10.00000
        assert outStr[::-1].find(".") == n_decimals

    def test_unit_rounding(self, real_usys):
        parameter = 0.001 * u.Unit("nm")
        n_decimals = 5
        outStr = real_usys.convert_parameter(parameter, n_decimals=n_decimals)
        assert outStr[::-1].find(".") == n_decimals

    def test_unitsystem_setup(self, real_usys):
        assert real_usys.usystem.name == "lammps_real"

        usys = LAMMPS_UnitSystems("lj", registry=u.UnitRegistry())
        assert usys.usystem.name == "lj"

    def test_dimensions_to_energy(self, real_usys):
        real_usys = LAMMPS_UnitSystems("real")
        parameter = 1 * u.kJ / u.nm * u.g
        # Note: parameter cannot divide out mass or time from energy
        out_parameter = real_usys._dimensions_to_energy(
            parameter.units.dimensions
        )
        assert str(out_parameter) == "(energy)*(mass)/(length)"

    def test_dimensions_to_charge(self, real_usys):
        parameter = 1 * u.coulomb / u.nm
        out_parameter = real_usys._dimensions_to_charge(
            parameter.units.dimensions
        )
        assert str(out_parameter) == "(charge)/(length)"

    def test_dimensions_thermal(self, real_usys):
        parameter = 1 * u.K
        out_parameter = real_usys._dimensions_from_thermal_to_energy(
            parameter.units.dimensions
        )
        assert str(out_parameter) == "(energy)"

    def test_get_dimensions(self):
        usys = LAMMPS_UnitSystems("electron")
        parametersList = list(
            map(
                lambda x: 1 * u.Unit(x, registry=usys.reg),
                [
                    "nm",
                    "kJ",
                    "kJ/mol",
                    "K",
                    "degree/angstrom",
                    "elementary_charge/mm",
                    "dimensionless",
                    "kg*m**2/s**2",
                    "coulomb",
                    "kcal/nm**2",
                    "K/nm",
                ],
            )
        )

        output_dimensionsList = [
            "length",
            "energy",
            "energy",
            "temperature",
            "angle/length",
            "charge/length",
            "dimensionless",
            "energy",
            "charge",
            "energy/length**2",
            "energy/length",
        ]
        for parameter, dim in zip(parametersList, output_dimensionsList):
            if str(parameter.units) == "K/nm":
                thermalize = True
            else:
                thermalize = False
            dimsStr = str(
                usys._get_output_dimensions(
                    parameter.units.dimensions, thermalize
                )
            )
            remove_parStr = dimsStr.translate({ord(i): None for i in "()"})
            assert remove_parStr == str(dim)

    def test_convert_parameters(self, typed_ethane, real_usys):
        parameter = typed_ethane.sites[0].atom_type.parameters["epsilon"]
        assert real_usys.convert_parameter(parameter) == "0.066"
        real_usys.usystem["energy"] = u.kJ / u.mol
        assert (
            real_usys.convert_parameter(parameter, n_decimals=6) == "0.276144"
        )
        usys = LAMMPS_UnitSystems("real")
        parameter = typed_ethane.bonds[0].bond_type.parameters["k"]
        assert usys.convert_parameter(parameter, n_decimals=0) == "680"
        parameter = typed_ethane.angles[0].angle_type.parameters["theta_eq"]
        assert usys.convert_parameter(parameter, name="theta_eq") == "110.700"
        parameter = typed_ethane.dihedrals[0].dihedral_type.parameters["c0"]
        assert usys.convert_parameter(parameter) == "0.150"

    def test_get_parameter_dimension(self):
        from gmso.utils.units import get_parameter_dimension

        assert (
            get_parameter_dimension(1 * u.kJ / u.mol / u.nm, "(energy)")
            == u.kJ / u.mol
        )
        assert get_parameter_dimension(1 * u.kJ / u.nm, "(length)") == u.nm
        assert (
            get_parameter_dimension(1 * u.kJ / u.mol / u.nm, "(length)") == u.nm
        )

    def test_convert_to_unit_system(self):
        # TODO: discuss if we would want a function to convert a whole
        # topology at once.
        # convert a whole topology to a specific unit system
        pass

    def test_generate_unit_styles(self):
        # TODO: write all unit systems for these engines.
        # look at libary of unit styles for lammps, gromacs, hoomd, gomc
        pass

    def test_lj_units(self, typed_ethane):
        # write out unit styles from ljUnitSystem and a dictonary of non-dimesnional values
        lj_usys = LAMMPS_UnitSystems("lj")
        bond_parameter = typed_ethane.bonds[0].bond_type.parameters["k"]
        errorStr = (
            "Missing conversion_factorDict for a dimensionless unit system."
        )
        with pytest.raises(ValueError, match=errorStr):
            lj_usys.convert_parameter(
                bond_parameter, conversion_factorDict=None
            )
        cfactorDict = {"energy": 0.276144 * u.kJ / u.mol, "length": 0.35 * u.nm}
        errorStr = f"Missing dimensionless constant in conversion_factorDict {cfactorDict}"
        with pytest.raises(ValueError, match=re.escape(errorStr)):
            lj_usys.convert_parameter(
                bond_parameter, conversion_factorDict=cfactorDict
            )
        cfactorDict["charge"] = 1 * u.coulomb
        cfactorDict["mass"] = 12.011 * u.amu
        outStr = lj_usys.convert_parameter(
            bond_parameter, conversion_factorDict=cfactorDict
        )
        assert outStr == str(round(284512.0 / 0.276144 * 0.35**2, 3))
        outStr = lj_usys.convert_parameter(
            typed_ethane.sites[0].atom_type.mass,
            conversion_factorDict=cfactorDict,
        )
        assert outStr == f"{(round(1.000, 3)):.3f}"

    def test_get_units(self, typed_ethane, real_usys):
        # get the unit system used for a topology
        typed_ethane.usystem = real_usys
        assert typed_ethane.usystem == real_usys

    def test_charmm_weighting_factors(self):
        # write out dihedrals while taking into account weighting
        pass

    def test_parameter_and_units_writing(self, real_usys):
        from gmso.utils.units import write_out_parameter_and_units

        x = 1 * u.kJ / u.mol
        outStr = write_out_parameter_and_units("x", x, real_usys)
        assert outStr == "x (kcal/mol)"

        x = 1 * u.rad
        outStr = write_out_parameter_and_units("theta_eq", x, real_usys)
        assert outStr == "theta_eq (degrees)"

        lj_usys = LAMMPS_UnitSystems("lj")
        outStr = write_out_parameter_and_units("x", x, lj_usys)
        assert outStr == "x (dimensionless)"
