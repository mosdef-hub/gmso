
import warnings
from typing import Optional, Union, List, Callable, Awaitable
from pydantic import ConfigDict, Field, field_serializer, field_validator
import unyt as u
from gmso.abc.abstract_site import Site
from gmso.core.atom import Atom
from gmso.core.virtual_type import VirtualPositionType, VirtualPotentialType

class VirtualSite(Site):
    """A generalized virtual site class in GMSO.

    Virtual sites are massless particles that represent off-atom charge sites, lone pairs, or other non-physical sites.
    
>>>>>>> 7d3a393 (added virtual site class attributes and added virtual type class)
    Attributes
    ----------
    charge : float
        The charge of the virtual site in elementary charge units.
    parent_atoms : List[Site]
        The real constituent atoms that define the virtual site's position.
    virtual_type:

    position:    
    """
    charge_: Optional[Union[u.unyt_quantity, float]] = Field(
        None, description="Charge of the virtual site", alias="charge"
    )
 
    parent_atoms_: List[Atom] = Field(
        ...,
        description="The 3 atoms involved in the angle.",
        alias="connection_members",
    ) 

    position_: Callable = Field(
       None, description="", alias="position"
   )


    virtual_type_: Optional[VirtualPositionType] = Field(
        default=None,
        description="virtual type for a virtual site .",
        alias="virtual_type",
    )
