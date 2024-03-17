"""Support non-bonded interactions between sites."""

import warnings
from typing import Optional, Set, Union

import unyt as u
from pydantic import ConfigDict, Field, field_serializer, field_validator

from gmso.abc.serialization_utils import unyt_to_dict
from gmso.core.parametric_potential import ParametricPotential
from gmso.utils._constants import UNIT_WARNING_STRING
from gmso.utils.expression import PotentialExpression
from gmso.utils.misc import (
    ensure_valid_dimensions,
    unyt_compare,
    unyt_to_hashable,
)
from gmso.utils.units import GMSO_UnitRegistry


class AtomType(ParametricPotential):
    """A description of non-bonded interactions between sites.

    This is a subclass of the gmso.core.Potential superclass.

    AtomType represents an atom type and includes the functional form
    describing its interactions and, optionally, other properties such as mass
    and charge.  This class inhereits from Potential, which stores the
    non-bonded interaction between atoms or sites. The functional form of the
    potential is stored as a `sympy` expression and the parameters, with units,
    are stored explicitly.
    """

    mass_: Optional[u.unyt_array] = Field(
        0.0 * u.gram / u.mol,
        description="The mass of the atom type",
        alias="mass",
    )

    charge_: Optional[u.unyt_array] = Field(
        0.0 * u.elementary_charge,
        description="The charge of the atom type",
        alias="charge",
    )

    atomclass_: Optional[str] = Field(
        "", description="The class of the atomtype", alias="atomclass"
    )

    doi_: Optional[str] = Field(
        "",
        description="Digital Object Identifier of publication where this atom type was introduced",
        alias="doi",
    )

    overrides_: Optional[Set[str]] = Field(
        set(),
        description="Set of other atom types that this atom type overrides",
        alias="overrides",
    )

    definition_: Optional[str] = Field(
        "",
        description="SMARTS string defining this atom type",
        alias="definition",
    )

    description_: Optional[str] = Field(
        "", description="Description for the AtomType", alias="description"
    )
    model_config = ConfigDict(
        alias_to_fields=dict(
            **ParametricPotential.model_config["alias_to_fields"],
            **{
                "mass": "mass_",
                "charge": "charge_",
                "atomclass": "atomclass_",
                "doi": "doi_",
                "overrides": "overrides_",
                "definition": "definition_",
                "description": "description_",
            },
        ),
    )

    def __init__(
        self,
        name="AtomType",
        mass=0.0 * u.gram / u.mol,
        charge=0.0 * u.elementary_charge,
        expression=None,
        parameters=None,
        potential_expression=None,
        independent_variables=None,
        atomclass="",
        doi="",
        overrides=None,
        definition="",
        description="",
        tags=None,
    ):
        if overrides is None:
            overrides = set()

        super(AtomType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables,
            potential_expression=potential_expression,
            mass=mass,
            charge=charge,
            atomclass=atomclass,
            doi=doi,
            overrides=overrides,
            description=description,
            definition=definition,
            tags=tags,
        )

    @property
    def charge(self):
        """Return the charge of the atom_type."""
        return self.__dict__.get("charge_")

    @property
    def mass(self):
        """Return the mass of the atom_type."""
        return self.__dict__.get("mass_")

    @property
    def atomclass(self):
        """Return the atomclass of the atom_type."""
        return self.__dict__.get("atomclass_")

    @property
    def doi(self):
        """Return the doi of the atom_type."""
        return self.__dict__.get("doi_")

    @property
    def overrides(self):
        """Return the overrides of the atom_type."""
        return self.__dict__.get("overrides_")

    @property
    def description(self):
        """Return the description of the atom_type."""
        return self.__dict__.get("description_")

    @property
    def definition(self):
        """Return the SMARTS string of the atom_type."""
        return self.__dict__.get("definition_")

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

    def clone(self, fast_copy=False):
        """Clone this AtomType, faster alternative to deepcopying."""
        return AtomType(
            name=str(self.name),
            tags=self.tags,
            expression=None,
            parameters=None,
            independent_variables=None,
            potential_expression=self.potential_expression.clone(fast_copy),
            mass=u.unyt_quantity(self.mass.value, self.mass.units),
            charge=u.unyt_quantity(self.charge.value, self.charge.units),
            atomclass=self.atomclass,
            doi=self.doi,
            overrides=(
                set(o for o in self.overrides) if self.overrides else None
            ),
            description=self.description,
            definition=self.definition,
        )

    def __hash__(self):
        """Return the unique hash of the object."""
        return id(self)

    def __eq__(self, other):
        if other is self:
            return True
        if not isinstance(other, AtomType):
            return False
        return (
            self.name == other.name
            and self.expression == other.expression
            and self.independent_variables == other.independent_variables
            and self.parameters.keys() == other.parameters.keys()
            and unyt_compare(
                self.parameters.values(), other.parameters.values()
            )
            and self.charge == other.charge
            and self.atomclass == other.atomclass
            and self.mass == other.mass
            and self.doi == other.doi
            and self.overrides == other.overrides
            and self.definition == other.definition
            and self.description == other.description
        )

    def _etree_attrib(self):
        attrib = super()._etree_attrib()
        if self.overrides == set():
            attrib.pop("overrides")
        mass = eval(attrib["mass"])
        charge = eval(attrib["charge"])

        attrib["mass"] = str(mass["array"])
        attrib["charge"] = str(charge["array"])

        return attrib

    def __repr__(self):
        """Return a formatted representation of the atom type."""
        desc = (
            f"<{self.__class__.__name__} {self.name},\n "
            f"expression: {self.expression},\n "
            f"id: {id(self)},\n "
            f"atomclass: {self.atomclass}>"
        )
        return desc

    @field_validator("mass_", mode="before")
    @classmethod
    def validate_mass(cls, mass):
        """Check to see that a mass is a unyt array of the right dimension."""
        default_mass_units = u.gram / u.mol
        if not isinstance(mass, u.unyt_array):
            warnings.warn(UNIT_WARNING_STRING.format("Masses", "g/mol"))
            mass *= u.gram / u.mol
        else:
            ensure_valid_dimensions(mass, default_mass_units)

        return mass

    @field_validator("charge_", mode="before")
    @classmethod
    def validate_charge(cls, charge):
        """Check to see that a charge is a unyt array of the right dimension."""
        if not isinstance(charge, u.unyt_array):
            warnings.warn(
                UNIT_WARNING_STRING.format("Charges", "elementary charge")
            )
            charge *= u.Unit(
                "elementary_charge", registry=GMSO_UnitRegistry().reg
            )
        else:
            ensure_valid_dimensions(charge, u.elementary_charge)

        return charge

    @staticmethod
    def _default_potential_expr():
        return PotentialExpression(
            expression="4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            independent_variables={"r"},
            parameters={
                "sigma": 0.3 * u.nm,
                "epsilon": 0.3 * u.Unit("kJ"),
            },
        )
