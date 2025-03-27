from typing import Callable, List, Optional, Union

import unyt as u
from pydantic import Field

from gmso.abc.abstract_site import Site
from gmso.core.virtual_type import VirtualSiteType


class VirtualSite(Site):
    """A generalized virtual site class in GMSO.

    Virtual sites are massless particles that represent off-atom charge sites, lone pairs, or other non-physical sites.

    Attributes
    ----------
    charge : float
        The charge of the virtual site in elementary charge units.
    parent_atoms : List[Site]
        The real constituent atoms that define the virtual site's position.
    virtual_type:

    position:
    """

    parent_atoms_: List[Site] = Field(
        ...,
        description="The parent atoms of the virtual site.",
        alias="parent_atoms",
    )

    charge_: Optional[Union[u.unyt_quantity, float]] = Field(
        None, description="Charge of the virtual site", alias="charge"
    )

    position_: Callable = Field(None, description="", alias="position")

    virtual_type_: Optional[VirtualSiteType] = Field(
        default=None,
        description="virtual type for a virtual site .",
        alias="virtual_type",
    )

    @property
    def parent_atoms(self) -> List[Site]:
        return self.__dict__.get("parent_atoms_")
