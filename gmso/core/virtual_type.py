import logging
from typing import Callable, Optional, Tuple, Union

import unyt as u
from pydantic import ConfigDict, Field, field_serializer, field_validator

from gmso.abc.gmso_base import GMSOBase
from gmso.abc.serialization_utils import unyt_to_dict
from gmso.core.parametric_potential import ParametricPotential
from gmso.utils._constants import UNIT_WARNING_STRING
from gmso.utils.expression import PotentialExpression
from gmso.utils.misc import ensure_valid_dimensions, unyt_compare
from gmso.utils.units import GMSO_UnitRegistry

logger = logging.getLogger(__name__)


class VirtualPositionType(ParametricPotential):
    """A description of the interaction between virtual partners.

    This is a subclass of the gmso.core.Potential superclass.

    VirtualPositionType represents a virtual type and includes the functional form
    describing the position of the virtual site from its parent atoms.
    The functional form of the potential is stored
    as a `sympy` expression and the parameters, with units, are stored
    explicitly. The AtomTypes that are used to define the virtual type are
    stored as `member_types`. The expression uses these atom_types with their
    positions as ri, rj, rk, ... , matching the number of atoms used to define
    the virtual site.

    Notes
    ----
    Inherits many functions from gmso.ParametricPotential:
        __eq__, _validate functions
    """

    model_config = ConfigDict(
        alias_to_fields=dict(
            **ParametricPotential.model_config["alias_to_fields"],
            **{},
        ),
    )

    def __init__(
        self,
        name="VirtualPositionType",
        expression=None,
        parameters=None,
        independent_variables=None,
        potential_expression=None,
        tags=None,
    ):
        super(VirtualPositionType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables,
            potential_expression=potential_expression,
            tags=tags,
        )

    @staticmethod
    def _default_potential_expr():
        return PotentialExpression(
            expression="ri + b*(rj-ri+a*(rk-rj))",
            parameters={
                "b": 1 * u.dimensionless,
                "a": 180 * u.dimensionless,
            },
            independent_variables={"ri", "rj", "rk"},
        )

    def clone(self, fast_copy=False):
        """Clone this VirtualPositionType, faster alternative to deepcopying."""
        return VirtualPositionType(
            name=str(self.name),
            tags=self.tags,
            expression=None,
            parameters=None,
            independent_variables=None,
            potential_expression=self.potential_expression.clone(fast_copy),
        )


class VirtualPotentialType(ParametricPotential):
    """A description of non-bonded interactions between sites.

    This is a subclass of the gmso.core.Potential superclass.

    VirtualPotentialType represents a virtual site type and includes the functional form
    describing its nonbonded interactions. This class inhereits from ParametricPotential, which stores the
    non-bonded interaction between atoms or sites. The functional form of the
    potential is stored as a `sympy` expression and the parameters, with units,
    are stored explicitly.
    """

    model_config = ConfigDict(
        alias_to_fields=dict(
            **ParametricPotential.model_config["alias_to_fields"],
            **{},
        ),
    )

    def __init__(
        self,
        name="VirtualPotentialType",
        expression=None,
        parameters=None,
        potential_expression=None,
        independent_variables=None,
        tags=None,
    ):
        super(VirtualPotentialType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables,
            potential_expression=potential_expression,
            tags=tags,
        )

    def clone(self, fast_copy=False):
        """Clone this VirtualPotentialType, faster alternative to deepcopying."""
        return VirtualPotentialType(
            name=str(self.name),
            tags=self.tags,
            expression=None,
            parameters=None,
            independent_variables=None,
            potential_expression=self.potential_expression.clone(fast_copy),
        )

    def __hash__(self):
        """Return the unique hash of the object."""
        return id(self)

    def __eq__(self, other):
        if other is self:
            return True
        if not isinstance(other, VirtualPotentialType):
            return False
        return (
            self.name == other.name
            and self.expression == other.expression
            and self.independent_variables == other.independent_variables
            and self.parameters.keys() == other.parameters.keys()
            and unyt_compare(self.parameters.values(), other.parameters.values())
        )

    def __repr__(self):
        """Return a formatted representation of the atom type."""
        desc = (
            f"<{self.__class__.__name__} {self.name},\n "
            f"expression: {self.expression},\n "
            f"id: {id(self)},\n "
        )
        return desc

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


class VirtualType(GMSOBase):
    """A generalized virtual site type class in GMSO.

    Virtual sites are massless particles that represent off-atom interaction/charge sites, lone pairs, or other non-physical sites.

    Attributes
    ----------
    name : str
        A generalized name for the type of the virtual site.
        See https://manual.gromacs.org/current/reference-manual/functions/interaction-methods.html#id3 for some Gromacs examples.
    member_types_: List[Str]
        The atom types of the constituent atoms that define the virtual site's position.
    member_classes_: List[Str]
        The atom classes of the constituent atoms that define the virtual site's position.
    charge : u.unyt_array
        The charge of the virtual site in elementary charge units.
    virtual_potential: gmso.core.virtual_type.VirtualPositionType
        The ParametricPotential that takes the specific `parent_atoms` and is used to generate the specific virtual site position.
    virtual_position: gmso.core.virtual_type.VirtualPotentialType
        The ParametricPotential that is used to store the nonbonded interactions of the virtual site type.
    """

    name_: str = Field(
        "VirtualType-3fd", description="Identifying Virtual Type", alias="name"
    )
    member_types_: Optional[Tuple[str, ...]] = Field(
        None,
        description="List-like of gmso.AtomType.name "
        "defining the members of this angle type",
        alias="member_types",
        min_length=0,
        max_length=12,
    )

    member_classes_: Optional[Tuple[str, ...]] = Field(
        None,
        description="List-like of gmso.AtomType.atomclass "
        "defining the members of this angle type",
        alias="member_classes",
        min_length=0,
        max_length=12,
    )
    identifier_: Optional[str] = Field(
        None,
        description="ForceField key identifier",
        alias="identifier",
    )

    charge_: Optional[u.unyt_array] = Field(
        0.0 * u.elementary_charge,
        description="The charge of the atom type",
        alias="charge",
    )

    position_: Callable = Field(None, description="", alias="position")

    virtual_position_: Optional[VirtualPositionType] = Field(
        default=VirtualPositionType(),
        description="virtual type for a virtual site.",
        alias="virtual_position",
    )

    virtual_potential_: Optional[VirtualPotentialType] = Field(
        default=VirtualPotentialType(),
        description="virtual type for a virtual site.",
        alias="virtual_potential",
    )

    doi_: Optional[str] = Field(
        "",
        description="Digital Object Identifier of publication where this atom type was introduced",
        alias="doi",
    )
    model_config = ConfigDict(
        alias_to_fields=dict(
            **ParametricPotential.model_config["alias_to_fields"],
            **{
                "charge": "charge_",
                "virtual_position": "virtual_position_",
                "virtual_potential": "virtual_potential_",
                "doi": "doi_",
                "member_classes": "member_classes_",
                "member_types": "member_types_",
                "identifier": "identifier_",
            },
        ),
    )

    @property
    def name(self) -> str:
        return self.__dict__.get("name_")

    @property
    def virtual_position(self) -> VirtualPositionType:
        return self.__dict__.get("virtual_position_")

    @property
    def virtual_potential(self) -> VirtualPotentialType:
        return self.__dict__.get("virtual_potential_")

    @property
    def member_types(self):
        return self.__dict__.get("member_types_")

    @property
    def member_classes(self):
        return self.__dict__.get("member_classes_")

    @property
    def identifier(self):
        return self.__dict__.get("identifier_")

    @property
    def charge(self):
        return self.__dict__.get("charge_")

    @property
    def doi(self):
        """Return the doi of the virtual_type."""
        return self.__dict__.get("doi_")

    @field_validator("charge_", mode="before")
    @classmethod
    def validate_charge(cls, charge):
        """Check to see that a charge is a unyt array of the right dimension."""
        if not isinstance(charge, u.unyt_array):
            logger.info(UNIT_WARNING_STRING.format("Charges", "elementary charge"))
            charge *= u.Unit("elementary_charge", registry=GMSO_UnitRegistry().reg)
        else:
            ensure_valid_dimensions(charge, u.elementary_charge)

        return charge

    @field_serializer("charge_")
    def serialize_charge(self, charge_: Union[u.unyt_quantity, None]):
        if charge_ is None:
            return None
        else:
            return unyt_to_dict(charge_)

    def __eq__(self, other):
        if other is self:
            return True
        if not isinstance(other, VirtualType):
            return False
        return (
            self.name == other.name
            and self.virtual_position == other.virtual_position
            and self.virtual_potential == other.virtual_potential
            and self.member_types == other.member_types
            and self.member_classes == other.member_classes
            and self.identifier == other.identifier
            and self.doi == other.doi
            and self.charge == other.charge
        )

    def clone(self, fast_copy=False):
        """Clone this VirtualType."""
        return VirtualType(
            name=str(self.name),
            virtual_position=self.virtual_position.clone(fast_copy),
            virtual_potential=self.virtual_potential.clone(fast_copy),
            member_types=self.member_types,
            member_classes=self.member_classes,
            identifier=self.identifier,
            charge=self.charge,
            doi=self.doi,
        )

    def get_parameters(self, copy=False):
        """Return parameters for this VirtualPotentialType and VirtualPositionType."""
        parametersDict = {}
        if copy:
            if self.virtual_potential:
                parametersDict["potential"] = {
                    k: u.unyt_quantity(v.value, v.units)
                    for k, v in self.virtual_potential.parameters.items()
                }
            if self.virtual_position:
                parametersDict["position"] = {
                    k: u.unyt_quantity(v.value, v.units)
                    for k, v in self.virtual_position.parameters.items()
                }
        else:
            if self.virtual_potential:
                parametersDict["potential"] = self.virtual_potential.parameters
            if self.virtual_position:
                parametersDict["position"] = self.virtual_position.parameters

        return parametersDict
