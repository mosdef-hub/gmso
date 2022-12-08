import warnings
from typing import Optional, Union

import unyt as u
from pydantic import Field, validator

from gmso.abc.abstract_site import Site
from gmso.core.atom_type import AtomType
from gmso.core.element import Element
from gmso.utils._constants import UNIT_WARNING_STRING
from gmso.utils.misc import ensure_valid_dimensions


class Bead(Site):
    __base_doc__ = """A bead represents a single coarse-grained association in a topology.

    Beads are the representation of an object within `gmso` that describes any general
    non-atomistic particle in gmso. Beads also contain information that are unique to
    elements vs other types of interaction sites in molecular simulations.

    Notes
    -----
    Beads have all the attributes inherited from the base Site class,
    The order of precedence when attaining properties `charge` and `mass` is:

        1. atom.charge > atom.atom_type.charge
        2. atom.mass > atom.atom_type.mass

    Examples
    --------
        >>> from gmso.core.bead import Bead
        >>> bead1 = Bead(name='SP1')

    See Also
    --------
    gmso.abc.AbstractSite
        An Abstract Base class for implementing site objects in GMSO. The class Bead bases from
        the gmso.abc.abstract site class
    """
    charge_: Optional[Union[u.unyt_quantity, float]] = Field(
        None,
        description="Charge of the bead",
    )

    mass_: Optional[Union[u.unyt_quantity, float]] = Field(
        None, description="Mass of the bead"
    )

    bead_type_: Optional[AtomType] = Field(
        None, description="AtomType associated with the atom"
    )

    @property
    def charge(self) -> Union[u.unyt_quantity, None]:
        charge = self.__dict__.get("charge_", None)
        atom_type = self.__dict__.get("atom_type_", None)
        if charge is not None:
            return charge
        elif bead_type is not None:
            return bead_type.charge
        else:
            return None

    @property
    def mass(self) -> Union[u.unyt_quantity, None]:
        mass = self.__dict__.get("mass_", None)
        atom_type = self.__dict__.get("atom_type_", None)
        if mass is not None:
            return mass
        elif bead_type is not None:
            return bead_type.mass
        else:
            return None

    @property
    def bead_type(self) -> Union[AtomType, None]:
        return self.__dict__.get("bead_type_", None)

    def __le__(self, other):
        if isinstance(other, Bead):
            return hash(self) <= hash(other)
        else:
            raise TypeError(
                f"Cannot compare equality between {type(self)} and {type(other)}"
            )

    def __lt__(self, other):
        if isinstance(other, Bead):
            return hash(self) < hash(other)
        else:
            raise TypeError(
                f"Cannot compare equality between {type(self)} and {type(other)}"
            )

    @validator("charge_")
    def is_valid_charge(cls, charge):
        if charge is None:
            return None
        if not isinstance(charge, u.unyt_array):
            warnings.warn(
                UNIT_WARNING_STRING.format("Charges", "elementary charge")
            )
            charge *= u.elementary_charge
        else:
            ensure_valid_dimensions(charge, u.elementary_charge)

        return charge

    @validator("mass_")
    def is_valid_mass(cls, mass):
        if mass is None:
            return None
        default_mass_units = u.gram / u.mol
        if not isinstance(mass, u.unyt_array):
            warnings.warn(UNIT_WARNING_STRING.format("Masses", "g/mol"))
            mass *= default_mass_units
        else:
            ensure_valid_dimensions(mass, default_mass_units)
        return mass

    class Config:
        extra = "forbid"

        fields = {
            "charge_": "charge",
            "mass_": "mass",
            "bead_type_": "bead_type",
        }

        alias_to_fields = {
            "charge": "charge_",
            "mass": "mass_",
            "bead_type": "bead_type_",
        }

        validate_assignment = True
