from typing import Callable, List, Optional, Union

import unyt as u
from pydantic import ConfigDict, Field

from gmso.abc.abstract_site import Site
from gmso.core.virtual_type import VirtualType
from gmso.exceptions import MissingPotentialError, NotYetImplementedWarning


class VirtualSite(Site):
    """A generalized virtual site class in GMSO.

    Virtual sites are massless particles that represent off-atom charge/interaction sites, lone pairs, or other non-physical sites.

    Attributes
    ----------
    charge : u.unyt_array
        The charge of the virtual site in elementary charge units. Will prioritize self.virtual_type.charge.
    parent_atoms : List[Site]
        The real constituent atoms that define the virtual site's position.
    virtual_type : gmso.core.virtual_type.VirtualType
        The type information, including parameters for virtual_position and virtual_potential, used to define
        the virtual site's interactions and positions
    """

    parent_atoms_: List[Site] = Field(
        default=[],
        description="The parent atoms of the virtual site.",
        alias="parent_atoms",
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
                "parent_atoms": "parent_atoms_",
            },
        ),
    )

    @property
    def parent_atoms(self) -> List[Site]:
        """Reminder that the order of atoms is fixed, such that atom index 1 corresponds to ri in the self.virtual_type.virtual_position expression."""
        return self.__dict__.get("parent_atoms_", [])

    def position(self):
        """Not yet implemented function to get position from virtual_type.virtual_position and parent_atoms."""
        if not self.virtual_potential:
            raise MissingPotentialError(
                "No VirtualType associated with this VirtualSite."
            )
        if not self.virtual_potential.virtual_position:
            raise MissingPotentialError(
                "No VirtualPositionType associated with this VirtualType."
            )
        # TODO: validate parent atoms matches virtual_type.virtual_position in terms of independent variables ri, rj, etc.
        # TODO: Generate position from atoms of parent_atoms and self.virtual_type.virtual_position.expression.
        raise NotYetImplementedWarning(
            "Need a functional to call from self.virtual_type.virtual_position, and plug in ri, rj, rk etc."
        )
        return 0.0

    def __repr__(self):
        return self.name + ": -".join(atom.__repr__() for atom in self.parent_atoms)

    @property
    def virtual_type(self):
        """Return the angle type if the angle is parametrized."""
        return self.__dict__.get("virtual_type_")
