import re

import pytest
import unyt as u

from gmso.tests.base_test import BaseTest
from gmso.utils.misc import unyt_to_hashable
from gmso.utils.units import LAMMPS_UnitSystems


class TestUnitHandling(BaseTest):
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

    def test_unit_conversion(self):
        usys = LAMMPS_UnitSystems("real")
        parameter = 0.001 * u.Unit("mm")
        n_decimals = 5
        outStr = usys.convert_parameter(parameter, n_decimals=n_decimals)
        assert float(outStr) == 10000.00000
        assert outStr[::-1].find(".") == n_decimals

        parameter = 1 * u.Unit("nm")
        n_decimals = 5
        outStr = usys.convert_parameter(parameter, n_decimals=n_decimals)
        assert float(outStr) == 10.00000
        assert outStr[::-1].find(".") == n_decimals

    def test_unit_rounding(self):
        usys = LAMMPS_UnitSystems("real")
        parameter = 0.001 * u.Unit("nm")
        n_decimals = 5
        outStr = usys.convert_parameter(parameter, n_decimals=n_decimals)
        assert outStr[::-1].find(".") == n_decimals

    def test_unitsystem_setup(self):
        usys = LAMMPS_UnitSystems("real")
        assert usys.system.name == "lammps_real"

        usys = LAMMPS_UnitSystems("lj", registry=u.UnitRegistry())
        assert usys.system.name == "lj"

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

    def test_convert_parameters(self, typed_ethane):
        parameter = typed_ethane.sites[0].atom_type.parameters["epsilon"]
        usys = LAMMPS_UnitSystems("real")
        assert usys.convert_parameter(parameter) == "0.066"
        usys.system["energy"] = u.kJ / u.mol
        assert usys.convert_parameter(parameter, n_decimals=6) == "0.276144"
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

    def test_overwrite_unit_output(self):
        # TODO: make a better usys.convert_parameter function to handle eq angles
        # write out the parameters that have specific outputs
        pass

    def test_generate_unit_styles(self):
        # TODO: write all unit systems for these engines.
        # look at libary of unit styles for lammps, gromacs, hoomd, gomc
        pass

    def test_lj_units(self, typed_ethane):
        # write out unit styles from ljUnitSystem and a dictonary of non-dimesnional values
        usys = LAMMPS_UnitSystems("lj")
        bond_parameter = typed_ethane.bonds[0].bond_type.parameters["k"]
        errorStr = (
            "Missing conversion_factorDict for a dimensionless unit system."
        )
        with pytest.raises(ValueError, match=errorStr):
            usys.convert_parameter(bond_parameter, conversion_factorDict=None)
        cfactorDict = {"energy": 0.276144 * u.kJ / u.mol, "length": 0.35 * u.nm}
        errorStr = f"Missing dimensionless constant in conversion_factorDict {cfactorDict}"
        with pytest.raises(ValueError, match=re.escape(errorStr)):
            usys.convert_parameter(
                bond_parameter, conversion_factorDict=cfactorDict
            )
        cfactorDict["charge"] = 1 * u.coulomb
        cfactorDict["mass"] = 12.011 * u.amu
        outStr = usys.convert_parameter(
            bond_parameter, conversion_factorDict=cfactorDict
        )
        assert outStr == str(round(284512.0 / 0.276144 * 0.35**2, 3))
        outStr = usys.convert_parameter(
            typed_ethane.sites[0].atom_type.mass,
            conversion_factorDict=cfactorDict,
        )
        assert outStr == f"{(round(1.000, 3)):.3f}"

    def test_get_units(self, typed_ethane):
        # get the unit system used for a topology
        usys = LAMMPS_UnitSystems("real")
        typed_ethane.system = usys
        assert typed_ethane.system == usys

    def test_charmm_weighting_factors(self):
        # write out dihedrals while taking into account weighting
        pass

    def test_parameter_and_units_writing(self):
        from gmso.utils.units import write_out_parameter_and_units

        base_unyts = LAMMPS_UnitSystems("real")
        x = 1 * u.kJ / u.mol
        outStr = write_out_parameter_and_units("x", x, base_unyts)
        assert outStr == "x (kcal/mol)"

        x = 1 * u.rad
        outStr = write_out_parameter_and_units("theta_eq", x, base_unyts)
        assert outStr == "theta_eq (degrees)"

        base_unyts = LAMMPS_UnitSystems("lj")
        outStr = write_out_parameter_and_units("x", x, base_unyts)
        assert outStr == "x (dimensionless)"
