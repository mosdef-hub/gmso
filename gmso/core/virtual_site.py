import string
from typing import Callable, List, Optional, Union

import unyt as u
from pydantic import ConfigDict, Field
from sympy import Matrix, symbols

from gmso.abc.abstract_site import Site, Molecule
from gmso.core.virtual_type import VirtualType
from gmso.exceptions import MissingPotentialError


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

    charge_: u.unyt_quantity = Field(
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

    def position(self) -> u.unyt_array:
        """On the fly position evaluation from virtual_type.virtual_position and parent_sites."""
        if not self.virtual_type:
            raise MissingPotentialError(
                "No VirtualType associated with this VirtualSite."
            )
        if not self.virtual_type.virtual_position:
            raise MissingPotentialError(
                "No VirtualPositionType associated with this VirtualType."
            )

        independent_namespace = {}
        for _, symbol in zip(range(len(self.parent_sites)), string.ascii_lowercase[8:]):
            x, y, z = symbols(f"r{symbol}1 r{symbol}2 r{symbol}3")
            independent_namespace[f"r{symbol}"] = Matrix([x, y, z])

        independent_parameters = {}
        for symbol, site in zip(string.ascii_lowercase[8:], self.parent_sites):
            for i, pos in enumerate(site.position):
                independent_parameters[f"r{symbol}{i + 1}"] = float(pos.value)

        # get units from parent sites
        unitsUnyt = self.parent_sites[0].position.units

        # perform expression evaluation
        return (
            self.virtual_type.virtual_position.potential_expression.evaluate(
                independent_namespace, independent_parameters
            )
            * unitsUnyt
        )

    def __repr__(self):
        return self.name + ": -".join(site.__repr__() for site in self.parent_sites)

    @property
    def virtual_type(self):
        """Return the virtual site type if the virtual site is parametrized."""
        return self.__dict__.get("virtual_type_")
    
    @property
    def charge(self) -> Union[u.unyt_quantity, None]:
        """Return the charge of the virtual site."""
        charge = self.__dict__.get("charge_", None)
        vtype = self.__dict__.get("virtual_type_", None)
        if charge is not None:
            return charge
        elif vtype is not None:
            return vtype.charge
        else:
            return 0.0 * u.elementary_charge

    @property
    def molecule(self) -> Union[Molecule, None]:
        """Return the molecule associated with the parent sites."""
        if not self.parent_sites:
            return None
        # Assume parent sites are part of the same molecule
        return self.parent_sites[0].molecule

