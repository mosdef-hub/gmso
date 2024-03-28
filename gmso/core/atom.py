"""Represent general atomic information in GMSO."""

import warnings
from typing import Optional, Union

import unyt as u
from pydantic import ConfigDict, Field, field_serializer, field_validator

from gmso.abc.abstract_site import Site
from gmso.abc.serialization_utils import unyt_to_dict
from gmso.core.atom_type import AtomType
from gmso.core.element import Element
from gmso.utils._constants import UNIT_WARNING_STRING
from gmso.utils.misc import ensure_valid_dimensions


class Atom(Site):
    """An atom represents a single element association in a topology.

    Atoms are the representation of an element within `gmso` that describes any general
    atom in a molecular simulation. Atoms also contain information that are unique to
    elements vs other types of interaction sites in molecular simulations.
    For example, charge, mass, and periodic table information.

    Notes
    -----
    Atoms have all the attributes inherited from the base Site class,
    The order of precedence when attaining properties `charge` and `mass` is:

    1. atom.charge > atom.atom_type.charge
    2. atom.mass > atom.atom_type.mass

    Examples
    --------
        >>> from gmso.core.atom import Atom
        >>> atom1 = Atom(name='lithium')

    See Also
    --------
    gmso.abc.AbstractSite
        An Abstract Base class for implementing site objects in GMSO. The class Atom bases from
        the gmso.abc.abstract site class
    """

    charge_: Optional[Union[u.unyt_quantity, float]] = Field(
        None, description="Charge of the atom", alias="charge"
    )

    mass_: Optional[Union[u.unyt_quantity, float]] = Field(
        None,
        description="Mass of the atom",
        alias="mass",
    )

    element_: Optional[Element] = Field(
        None,
        description="Element associated with the atom",
        alias="element",
    )

    atom_type_: Optional[AtomType] = Field(
        None, description="AtomType associated with the atom", alias="atom_type"
    )

    restraint_: Optional[dict] = Field(
        default=None,
        description="""
        Restraint for this atom, must be a dict with the following keys:
        'kx', 'ky', 'kz' (unit of energy/(mol * length**2),
        Refer to https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html
        for more information.
        """,
        alias="restraint",
    )

    model_config = ConfigDict(
        alias_to_fields=dict(
            **Site.model_config["alias_to_fields"],
            **{
                "charge": "charge_",
                "mass": "mass_",
                "element": "element_",
                "atom_type": "atom_type_",
                "restraint": "restraint_",
            },
        ),
    )

    @property
    def charge(self) -> Union[u.unyt_quantity, None]:
        """Return the charge of the atom."""
        charge = self.__dict__.get("charge_", None)
        atom_type = self.__dict__.get("atom_type_", None)
        if charge is not None:
            return charge
        elif atom_type is not None:
            return atom_type.charge
        else:
            return None

    @property
    def mass(self) -> Union[u.unyt_quantity, None]:
        """Return the mass of the atom."""
        mass = self.__dict__.get("mass_", property)
        atom_type = self.__dict__.get("atom_type_", None)
        if mass is not None:
            return mass
        elif atom_type is not None:
            return atom_type.mass
        else:
            return None

    @property
    def element(self) -> Union[Element, None]:
        """Return the element associated with the atom."""
        return self.__dict__.get("element_", None)

    @property
    def atom_type(self) -> Union[AtomType, property]:
        """Return the atom_type associated with the atom."""
        return self.__dict__.get("atom_type_", None)

    @property
    def restraint(self):
        """Return the restraint of this atom."""
        return self.__dict__.get("restraint_")

    @field_serializer("charge_")
    def serialize_charge(self, charge_: Union[u.unyt_quantity, None]):
        if charge_ is None:
            return None
        else:
            return unyt_to_dict(charge_)

    @field_serializer("mass_")
    def serialize_mass(self, mass_: Union[u.unyt_quantity, None]):
        if mass_ is None:
            return None
        else:
            return unyt_to_dict(mass_)

    @field_serializer("restraint_")
    def serialize_restraint(self, restraint_: Union[dict, None]):
        if restraint_ is None:
            return None
        else:
            converted_restraint = {
                key: unyt_to_dict(val) for key, val in restraint_.items()
            }
        return converted_restraint

    def clone(self):
        """Clone this atom."""
        return Atom(
            name=self.name,
            label=self.label,
            group=self.group,
            molecule=self.molecule,
            residue=self.residue,
            position=self.position,
            charge=self.charge,
            mass=self.mass,
            element=self.element,
            atom_type=(
                property if not self.atom_type else self.atom_type.clone()
            ),
        )

    def __le__(self, other):
        """Less than or equal to operator."""
        if isinstance(other, Atom):
            return hash(self) <= hash(other)
        else:
            raise TypeError(
                f"Cannot compare equality between {type(self)} and {type(other)}"
            )

    def __lt__(self, other):
        """Less than operator."""
        if isinstance(other, Atom):
            return hash(self) < hash(other)
        else:
            raise TypeError(
                f"Cannot compare equality between {type(self)} and {type(other)}"
            )

    @field_validator("charge_")
    @classmethod
    def is_valid_charge(cls, charge):
        """Ensure that the charge is physically meaningful."""
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

    @field_validator("mass_")
    @classmethod
    def is_valid_mass(cls, mass):
        """Ensure that the mass is physically meaningful."""
        if mass is None:
            return None
        default_mass_units = u.gram / u.mol
        if not isinstance(mass, u.unyt_array):
            warnings.warn(UNIT_WARNING_STRING.format("Masses", "g/mol"))
            mass *= default_mass_units
        else:
            ensure_valid_dimensions(mass, default_mass_units)
        return mass
