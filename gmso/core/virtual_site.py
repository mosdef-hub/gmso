from typing import Callable, List, Optional, Union

import unyt as u
from pydantic import ConfigDict, Field

from gmso.abc.abstract_site import Site
from gmso.core.virtual_type import VirtualType
from gmso.exceptions import GMSOError, MissingPotentialError


class VirtualSite(Site):
    """A generalized virtual site class in GMSO.

    Virtual sites are massless particles that represent off-atom charge/interaction sites, lone pairs, or other non-physical sites.

    Attributes
    ----------
    charge : u.unyt_array
        The charge of the virtual site in elementary charge units. Will prioritize self.virtual_type.charge.
    parent_sites : List[Site]
        The real constituent sites that define the virtual site's position.
    virtual_type : gmso.core.virtual_type.VirtualType
        The type information, including parameters for virtual_position and virtual_potential, used to define
        the virtual site's interactions and positions
    """

    parent_sites_: List[Site] = Field(
        default=[],
        description="The parent sites of the virtual site.",
        alias="parent_sites",
    )

    charge_: Optional[Union[u.unyt_quantity, float]] = Field(
        None, description="Charge of the virtual site", alias="charge"
    )

    position_: Callable = Field(None, description="", alias="position")

    virtual_type_: Optional[VirtualType] = Field(
        default=None,
        description="virtual type for a virtual site.",
        alias="virtual_type",
    )

    model_config = ConfigDict(
        alias_to_fields=dict(
            **Site.model_config["alias_to_fields"],
            **{
                "charge": "charge_",
                "virtual_type": "virtual_type_",
                "parent_sites": "parent_sites_",
            },
        ),
    )

    @property
    def parent_sites(self) -> List[Site]:
        """Reminder that the order of sites is fixed, such that site index 1 corresponds to ri in the self.virtual_type.virtual_position expression."""
        return self.__dict__.get("parent_sites_", [])

    def position(self) -> str:
        """Not yet implemented function to get position from virtual_type.virtual_position and parent_sites."""
        if not self.virtual_type:
            raise MissingPotentialError(
                "No VirtualType associated with this VirtualSite."
            )
        if not self.virtual_type.virtual_position:
            raise MissingPotentialError(
                "No VirtualPositionType associated with this VirtualType."
            )
        import string

        import numpy as np
        import sympy
        from sympy import Symbol, sympify

        class norm(sympy.Function):
            @classmethod
            def eval(cls, arg):
                return None

        # Evaluate vector norm
        def norm_evaluation(matrix_arg):
            return np.linalg.norm(matrix_arg)

        # String expression
        expr_string = str(self.virtual_type.virtual_position.expression)
        namespace = {
            sym: Symbol(sym) for sym in self.virtual_type.virtual_position.parameters
        }
        namespace["norm"] = norm
        args = [
            namespace[sym] for sym in self.virtual_type.virtual_position.parameters
        ]  # args to lambdify

        for atom, symbol in zip(
            range(len(self.parent_sites)), string.ascii_lowercase[8:]
        ):
            x, y, z = sympy.symbols(f"r{symbol}1 r{symbol}2 r{symbol}3")
            namespace[f"r{symbol}"] = sympy.Matrix([x, y, z])
            args.extend([x, y, z])

        # Parse expression
        try:
            expr = sympify(expr_string, locals=namespace)
        except (ValueError, TypeError):
            raise GMSOError(
                f"Expression passed {expr_string=} was not viable in sympy."
            )

        f = sympy.lambdify(args, expr, modules=[{"norm": norm_evaluation}, "numpy"])

        # Evaluate
        parameters = {
            param: val.to_value()
            for param, val in self.virtual_type.virtual_position.parameters.items()
        }
        for symbol, site in zip(string.ascii_lowercase[8:], self.parent_sites):
            for i, pos in enumerate(site.position):
                parameters[f"r{symbol}{i + 1}"] = float(pos.value)
        result = f(**parameters)
        if isinstance(result, Symbol) or not np.issubdtype(result.dtype, np.floating):
            raise GMSOError(
                f"{self=} expression was not able to be fully evaluated. Unknown parameters left defined in {expr_string=}. Position evaluation is: {result=}"
            )
        elif any(np.isnan(result)):
            raise GMSOError(
                f"Failed evaluation of {self=} with {expr_string=} and {parameters=}."
            )
        return np.array(result).T * u.nm

    def __repr__(self):
        return self.name + ": -".join(site.__repr__() for site in self.parent_sites)

    @property
    def virtual_type(self):
        """Return the virtual site type if the virtual site is parametrized."""
        return self.__dict__.get("virtual_type_")
