"""Source of available units registered within GMSO."""

import numpy as np
import unyt as u


class GMSO_UnitRegsitry(object):
    """A default unit registry class.

    The basic units that need to be added for various unit conversions done
    throughout GMSO.

    Attributes
    ----------
    reg : u.UnitRegistry
        The unit registry useful for conversions commonly used in molecular topologies
    """

    def __init__(self):
        self._reg = u.UnitRegistry()
        conversion = (
            1 * getattr(u.physical_constants, "elementary_charge").value
        )
        self.register_unit(
            self._reg[u.dimensions.current_mks, u.dimensions.time],
            conversion,
            "elementary_charge",
            r"\rm{e}",
        )

    @staticmethod
    def register_unit(
        registry,
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
        registry.add(
            "elementary_charge",
            base_value=conversion,
            dimensions=dim,
            tex_repr=tex_repr,
        )

    @classmethod
    def default_reg(cls):
        """Return a default registry with extra units defined outside of unyt.

        Returns
        -------
        reg : u.unyt_registy
            A unyt registry with commonly used conversions defined.
        """
        reg = u.UnitRegistry()
        conversion = (
            1 * getattr(u.physical_constants, "elementary_charge").value
        )
        dimensionsList = [u.dimensions.current_mks, u.dimensions.time]
        name = "elementary_charge"
        symbol = r"\rm{e}"
        cls.register_unit(reg, name, conversion, dimensionsList, symbol)
        return reg
