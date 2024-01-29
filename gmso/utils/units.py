"""Source of available units registered within GMSO."""

import re

import numpy as np
import unyt as u
from sympy import Symbol

from gmso.exceptions import NotYetImplementedWarning


class GMSO_UnitRegistry(object):
    """A default unit registry class.

    The basic units that need to be added for various unit conversions done
    throughout GMSO.

    Attributes
    ----------
    reg : u.UnitRegistry
        The unit registry useful for conversions commonly used in molecular topologies
    """

    def __init__(self):
        self.reg_ = u.UnitRegistry()
        register_general_units(self.reg)

    def register_unit(
        self,
        name: str,
        conversion: float,
        dimensionsList: list,
        tex_repr=None,
    ):
        """Add units to the self.reg UnitRegistry.

        Parameters
        ----------
        registry : u.unyt_registy, required
            Unit registry to add the unit to. See unyt.unyt_registry for more information
        dimensionsList : list, required
            A list of the dimensions that the unit will be registered under. If using the inverse of a dimension
            be sure to supply 1/u.dimension as the element of the list.
        conversion : float, required
            The numerical value for the conversion in SI units with the same dimensions. See unyt.unyt_registry.add
            module for more information
        name : str, required
            Then name of the unyt to be referenced as string when calling u.Unit("unit_name")
        tex_repr : str, optional, default None
            The latex representation that is used to visualze the unit when pretty print is used.


        """
        dim = np.prod(dimensionsList)
        if not tex_repr:
            tex_repr = r"\rm{name}"
        self.reg_.add(
            symbol=name,
            base_value=conversion,
            dimensions=dim,
            tex_repr=tex_repr,
        )

    @property
    def reg(self):
        """Return the UnitRegistry attribute for the class."""
        return self.__dict__.get("reg_")

    @staticmethod
    def default_reg():
        """Return a default registry with extra units defined outside of unyt.

        Returns
        -------
        reg : u.unyt_registy
            A unyt registry with commonly used conversions defined.
        """
        reg = u.UnitRegistry()
        register_general_units(reg)
        return reg


def register_general_units(reg: u.UnitRegistry):
    """Register units that are generally useful to a basic unyt UnitSystem."""
    elementary_charge_conversion = (
        1 * getattr(u.physical_constants, "elementary_charge").value
    )
    charge_dim = u.dimensions.current_mks * u.dimensions.time
    reg.add(
        "elementary_charge",
        elementary_charge_conversion,
        charge_dim,
        r"\rm{e}",
    )  # proton charge

    kb_conversion = (
        1 * getattr(u.physical_constants, "boltzmann_constant_mks").value
    )
    kb_dim = u.dimensions.energy / u.dimensions.temperature
    reg.add(
        "kb", base_value=kb_conversion, dimensions=kb_dim, tex_repr=r"\rm{kb}"
    )  # boltzmann temperature

    bohr_rad_conversion = (
        4
        * np.pi
        * getattr(u.physical_constants, "reduced_planck_constant").value ** 2
        * getattr(u.physical_constants, "eps_0").value
        / (
            getattr(u.physical_constants, "electron_charge").value ** 2
            * getattr(u.physical_constants, "electron_mass").value
        )
    )
    bohr_rad_dim = u.dimensions.length
    reg.add(
        "a0",
        base_value=bohr_rad_conversion,
        dimensions=bohr_rad_dim,
        tex_repr=r"\rm{a0}",
    )  # bohr radius

    hartree_conversion = (
        getattr(u.physical_constants, "reduced_planck_constant").value ** 2
        / u.Unit("a0", registry=reg).base_value ** 2
        / getattr(u.physical_constants, "electron_mass").value
    )
    hartree_dim = u.dimensions.energy
    reg.add(
        "Ehartree",
        base_value=hartree_conversion,
        dimensions=hartree_dim,
        tex_repr=r"\rm{Ehartree}",
    )  # Hartree energy

    static_coulomb_conversion = np.sqrt(
        10**9 / (4 * np.pi * getattr(u.physical_constants, "eps_0").value)
    )
    charge_dim = u.dimensions.charge
    reg.add(
        "Statcoulomb_charge",
        base_value=static_coulomb_conversion,
        dimensions=charge_dim,
        tex_repr=r"\rm{Statcoulomb_charge}",
    )  # Static charge


class LAMMPS_UnitSystems:
    """Set of a unit systems distributed in LAMMPS (https://docs.lammps.org/units.html)."""

    def __init__(self, style: str, registry=None):
        if registry:
            self.reg_ = registry
        else:
            self.reg_ = GMSO_UnitRegistry().reg
        self.usystem_ = self.usystem_from_str(styleStr=style, reg=self.reg_)

    @property
    def usystem(self):
        """Return the UnitSystem attribute for the class."""
        return self.__dict__.get("usystem_")

    @property
    def reg(self):
        """Return the UnytRegistry attribute for the class."""
        return self.__dict__.get("reg_")

    def usystem_from_str(self, styleStr: str, reg: u.UnitRegistry):
        """Get systems for unit style."""
        #  NOTE: the when an angle is measured in lammps is not straightforwards. It depends not on the unit_style, but on the
        # angle_style, dihedral_style, or improper_style. For examples, harmonic angles, k is specificed in energy/radian, but the
        # theta_eq is written in degrees. For fourier dihedrals, d_eq is specified in degrees. When adding new styles, make sure that
        # this behavior is accounted for when converting the specific potential_type in the function
        # _parameter_converted_to_float
        if styleStr == "real":
            base_units = u.UnitSystem(
                "lammps_real",
                length_unit="angstrom",
                mass_unit="amu",
                time_unit="fs",
                temperature_unit="K",
                angle_unit="rad",
                registry=reg,
            )
            base_units["energy"] = "kcal/mol"
            base_units["charge"] = "elementary_charge"
        elif styleStr == "metal":
            base_units = u.UnitSystem(
                "lammps_metal",
                length_unit="angstrom",
                mass_unit="amu",
                time_unit="picosecond",
                temperature_unit="K",
                angle_unit="rad",
                registry=reg,
            )
            base_units["energy"] = "eV"
            base_units["charge"] = "elementary_charge"
        elif styleStr == "si":
            base_units = u.UnitSystem(
                "lammps_si",
                length_unit="m",
                mass_unit="kg",
                time_unit="s",
                temperature_unit="K",
                angle_unit="rad",
                registry=reg,
            )
            base_units["energy"] = "joule"
            base_units["charge"] = "coulomb"
        elif styleStr == "cgs":
            base_units = u.UnitSystem(
                "lammps_cgs",
                length_unit="cm",
                mass_unit="g",
                time_unit="s",
                temperature_unit="K",
                angle_unit="rad",
                registry=reg,
            )
            base_units["energy"] = "erg"
            # Statcoulomb is strange. It is not a 1:1 correspondance to charge, with base units of
            # mass**1/2*length**3/2*time**-1.
            # However, assuming it is referring to a static charge and not a flux, it can be
            # converted to coulomb units. See the registry for the unit conversion to Coulombs
            base_units["charge"] = "Statcoulomb_charge"
        elif styleStr == "electron":
            base_units = u.UnitSystem(
                "lammps_electron",
                length_unit="a0",
                mass_unit="amu",
                time_unit="s",
                temperature_unit="K",
                angle_unit="rad",
                registry=reg,
            )
            base_units["energy"] = "Ehartree"
            base_units["charge"] = "elementary_charge"
        elif styleStr == "micro":
            base_units = u.UnitSystem(
                "lammps_micro",
                length_unit="um",
                mass_unit="picogram",
                time_unit="us",
                temperature_unit="K",
                angle_unit="rad",
                registry=reg,
            )
            base_units["energy"] = "ug*um**2/us**2"
            base_units["charge"] = "picocoulomb"
        elif styleStr == "nano":
            base_units = u.UnitSystem(
                "lammps_nano",
                length_unit="nm",
                mass_unit="attogram",
                time_unit="ns",
                temperature_unit="K",
                angle_unit="rad",
                registry=reg,
            )
            base_units["energy"] = "attogram*nm**2/ns**2"
            base_units["charge"] = "elementary_charge"
        elif styleStr == "lj":
            base_units = ljUnitSystem(reg)
        else:
            raise NotYetImplementedWarning

        return base_units

    def convert_parameter(
        self,
        parameter,
        conversion_factorDict=None,
        n_decimals=3,
        name="",
    ):
        """Take a given parameter, and return a string of the parameter in the given style.

        This function will check the base_unyts, which is a unyt.UnitSystem object,
        and convert the parameter to those units based on its dimensions. It can
        also generate dimensionless units via normalization from conversion_factorsDict.

        Parameters
        ----------
        parameter : unyt.array.unyt_quantity
            Any parameter to convert to a string in the dimensions of self.usystem
        conversion_factorDict : dict, default=None
            If the self.usystem is ljUnitSystem, handle conversion
        n_decimals : int, default=3
            The number of decimals used in string .f formatting
        name : string, default=""
            Additionally information about the parameter, required if handling a specific parameter
            differently than the default self.usystem would.

        Returns
        -------
        outStr : str
            The parameter converted via self.usystem, and foramtted as a float string.
        """
        if name in [
            "theta_eq",
            "chieq",
            "phi_eq",
        ]:  # eq angle are always in degrees
            return f"{round(float(parameter.to('degree').value), n_decimals):.{n_decimals}f}"
        new_dims = self._get_output_dimensions(parameter.units.dimensions)
        if isinstance(self.usystem, ljUnitSystem):
            if not conversion_factorDict:
                raise ValueError(
                    "Missing conversion_factorDict for a dimensionless unit system."
                )
            elif not np.all(
                [
                    key in conversion_factorDict
                    for key in ["energy", "length", "mass", "charge"]
                ]
            ):
                raise ValueError(
                    f"Missing dimensionless constant in conversion_factorDict {conversion_factorDict}"
                )
            # multiply object -> split into length, mass, energy, charge -> grab conversion factor from dict
            # first replace energy for (length)**2*(mass)/(time)**2 u.dimensions.energy. Then iterate through the free symbols
            # and figure out a way how to add those to the overall conversion factor
            dim_info = new_dims.as_terms()
            conversion_factor = 1
            for exponent, ind_dim in zip(dim_info[0][0][1][1], dim_info[1]):
                factor = conversion_factorDict.get(
                    ind_dim.name[1:-1],
                    1 * self.usystem[ind_dim.name[1:-1]],  # default value of 1
                )  # replace () in name
                current_unit = get_parameter_dimension(parameter, ind_dim.name)
                factor = factor.to(
                    current_unit
                )  # convert factor to units of parameter
                conversion_factor *= float(factor) ** (exponent)
            return f"""{round(
                float(parameter / conversion_factor),
                n_decimals
            ):.{n_decimals}f}"""  # Assuming that conversion factor is in right units
        new_dimStr = str(new_dims)
        ind_units = re.sub("[^a-zA-Z]+", " ", new_dimStr).split()
        for unit in ind_units:
            new_dimStr = new_dimStr.replace(unit, str(self.usystem[unit]))
        outFloat = float(
            parameter.to(u.Unit(new_dimStr, registry=self.usystem.registry))
        )

        return f"{outFloat:.{n_decimals}f}"

    @staticmethod
    def _dimensions_to_energy(dims):
        """Take a set of dimensions and substitute in Symbol("energy") where possible."""
        symsStr = str(dims.free_symbols)
        energy_inBool = np.all(
            [dimStr in symsStr for dimStr in ["time", "mass"]]
        )  # TODO: this logic could be improved, it might fail on complex
        # units where the units are energy/mass/time**2, or something where the
        # dimensions are cancelled out
        if not energy_inBool:
            return dims
        energySym = Symbol(
            "(energy)"
        )  # create dummy symbol to replace in equation
        dim_info = dims.as_terms()
        time_idx = np.where(
            list(map(lambda x: x.name == "(time)", dim_info[1]))
        )[0][0]
        energy_exp = (
            dim_info[0][0][1][1][time_idx] // 2
        )  # energy has 1/time**2 in it, so this is the hint of how many
        return (
            dims
            * u.dimensions.energy**energy_exp
            * energySym ** (-1 * energy_exp)
        )

    @staticmethod
    def _dimensions_to_charge(dims):
        """Take a set of dimensions and substitute in Symbol("charge") where possible."""
        symsStr = str(dims.free_symbols)
        charge_inBool = np.all(
            [dimStr in symsStr for dimStr in ["current_mks"]]
        )
        if not charge_inBool:
            return dims
        chargeSym = Symbol(
            "(charge)"
        )  # create dummy symbol to replace in equation
        dim_info = dims.as_terms()
        current_idx = np.where(
            list(map(lambda x: x.name == "(current_mks)", dim_info[1]))
        )[0][0]
        charge_exp = dim_info[0][0][1][1][
            current_idx
        ]  # charge has (current_mks) in it, so this is the hint of how many
        return (
            dims
            * u.dimensions.charge ** (-1 * charge_exp)
            * chargeSym**charge_exp
        )

    @staticmethod
    def _dimensions_from_thermal_to_energy(dims):
        """Take a set of dimensions and substitute in Symbol("energy") to replace temperature."""
        symsStr = str(dims.free_symbols)
        temp_inBool = np.all([dimStr in symsStr for dimStr in ["temperature"]])
        if not temp_inBool:
            return dims
        energySym = Symbol(
            "(energy)"
        )  # create dummy symbol to replace in equation
        dim_info = dims.as_terms()
        temp_idx = np.where(
            list(map(lambda x: x.name == "(temperature)", dim_info[1]))
        )[0][0]
        temp_exp = dim_info[0][0][1][1][
            temp_idx
        ]  # energy has 1/time**2 in it, so this is the hint of how many
        return (
            dims / u.dimensions.temperature**temp_exp * energySym ** (temp_exp)
        )

    @classmethod
    def _get_output_dimensions(cls, dims, thermal_equivalence=False):
        if str(dims) == "1":  # use string as all dims can be converted
            return u.dimensionless
        dims = cls._dimensions_to_energy(dims)
        dims = cls._dimensions_to_charge(dims)
        if thermal_equivalence:
            dims = cls._dimensions_from_thermal_to_energy(dims)
        return dims


class ljUnitSystem:
    """Use this so the empty unitsystem has getitem magic method."""

    def __init__(self, reg: u.UnitRegistry):
        self.registry = reg
        self.name = "lj"

    def __getitem__(self, item):
        """Return dimensionless units unless angle."""
        if item == "angle":
            return u.Unit("degree")
        return u.Unit("dimensionless")


def get_parameter_dimension(parameter, dimension):
    """Return a unit from the parameter in a given dimension."""
    param_terms = parameter.units.expr.as_terms()
    uStr = ""
    for symbol, exp in zip(param_terms[-1], param_terms[0][0][1][1]):
        outputDim = LAMMPS_UnitSystems._get_output_dimensions(
            u.Unit(symbol).dimensions
        )
        if str(outputDim) == dimension:
            uStr += f"{symbol}*"
        elif (
            str(outputDim) == "dimensionless" and dimension == "(energy)"
        ):  # add mol to units of energy
            uStr += f"{symbol}**{exp}*"
        elif (
            str(outputDim) == "dimensionless" and dimension == "(mass)"
        ):  # add mol to mass amu
            uStr += f"{symbol}**{exp}*"
    return u.Unit(uStr[:-1])


def write_out_parameter_and_units(parameter_name, parameter, base_unyts=None):
    """Take a parameter and return a heading for the parameter and the units it should be in.

    Parameters
    ----------
    parameter_name : str
        The name of the unyt parameter to  be written. The dict key of the
        parameter associated with the GMSO object.
    parameter : unyt.array.unyt_quantity
        The unyt object with the units to be written out.
    base_unyts : LAMMPS_UnitSystem or a more general GMSO_UnitSystem
        The base units that house the relevant registry for
        converting parameters into a specified system.


    Returns
    -------
    output_parameter_units : str
        parameter with name converted into output unyt system. Useful for
        labeling parameters in output files, such as .data or .top files.
    """
    if not base_unyts:
        return f"{parameter_name} ({parameter.units})"
    if parameter_name in ["theta_eq", "phi_eq"]:
        return f"{parameter_name} ({'degrees'})"  # default to always degrees
    if base_unyts.usystem.name == "lj":
        return f"{parameter_name} ({'dimensionless'})"
    new_dims = LAMMPS_UnitSystems._get_output_dimensions(
        parameter.units.dimensions
    )
    new_dimStr = str(new_dims)
    ind_units = re.sub("[^a-zA-Z]+", " ", new_dimStr).split()
    for unit in ind_units:
        new_dimStr = new_dimStr.replace(unit, str(base_unyts.usystem[unit]))

    outputUnyt = str(
        parameter.to(u.Unit(new_dimStr, registry=base_unyts.reg)).units
    )
    return f"{parameter_name} ({outputUnyt})"


def convert_params_units(
    potentials,
    expected_units_dim,
    base_units,
):
    """Convert parameters' units in the potential to that specified in the base_units.

    Parameters
    ----------
    potentials : list
        Set of potentials to apply the conversion to.
    expected_units_dim : dict
        The dimensionality expected for all parameters in the potential. This allows
        the given dimensions to be converted to via the base_units specified in the
        unit system.
    base_units : dict
        The units to use for conversion. Must have keys of "length", "energy",
        and "mass". These are the base units any parameter will be converted into.

    Returns
    -------
    converted_potentials : list
        the input potentials converted into the base units given by
        base_units `dict`.
    """
    converted_potentials = list()
    for potential in potentials:
        converted_params = dict()
        for parameter in potential.parameters:
            unit_dim = expected_units_dim[parameter]
            ind_units = re.sub("[^a-zA-Z]+", " ", unit_dim).split()
            for unit in ind_units:
                unit_dim = unit_dim.replace(
                    unit,
                    f"({str(base_units[unit].value)} * {str(base_units[unit].units)})",
                )

            converted_params[parameter] = potential.parameters[parameter].to(
                unit_dim
            )
        potential.parameters = converted_params
        converted_potentials.append(potential)
    return converted_potentials
