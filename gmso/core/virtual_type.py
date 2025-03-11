
from typing import Optional, Tuple, Union, Set
import warnings
import unyt as u
from pydantic import ConfigDict, Field ,field_serializer, field_validator

from gmso.abc.serialization_utils import unyt_to_dict
from gmso.core.parametric_potential import ParametricPotential
from gmso.utils.expression import PotentialExpression
from gmso.utils._constants import UNIT_WARNING_STRING
from gmso.utils.misc import ensure_valid_dimensions, unyt_compare
from gmso.utils.units import GMSO_UnitRegistry

class VirtualPositionType(ParametricPotential):
    """A descripton of the interaction between 3 bonded partners.

    This is a subclass of the gmso.core.Potential superclass.

    AngleType represents an angle type and includes the functional form
    describing its interactions. The functional form of the potential is stored
    as a `sympy` expression and the parameters, with units, are stored
    explicitly.  The AtomTypes that are used to define the angle type are
    stored as `member_types`.

    Notes
    ----
    Inherits many functions from gmso.ParametricPotential:
        __eq__, _validate functions
    """

    member_types_: Optional[Tuple[str, str, str]] = Field(
        None,
        description="List-like of gmso.AtomType.name "
        "defining the members of this angle type",
        alias="member_types",
    )

    member_classes_: Optional[Tuple[str, str, str]] = Field(
        None,
        description="List-like of gmso.AtomType.atomclass "
        "defining the members of this angle type",
        alias="member_classes",
    )
    model_config = ConfigDict(
        alias_to_fields=dict(
            **ParametricPotential.model_config["alias_to_fields"],
            **{
                "member_types": "member_types_",
                "member_classes": "member_classes_",
            },
        ),
    )

    def __init__(
        self,
        name="VirtualPositionType",
        expression=None,
        parameters=None,
        independent_variables=None,
        potential_expression=None,
        member_types=None,
        member_classes=None,
        tags=None,
    ):
        super(VirtualPositionType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables,
            potential_expression=potential_expression,
            member_types=member_types,
            member_classes=member_classes,
            tags=tags,
        )

    @staticmethod
    def _default_potential_expr():
        return PotentialExpression(
            expression="ri + b*(rj-ri+a*(rk-rj))/(mag rijk)",
            parameters={
                "k": 1000 * u.Unit("kJ / (deg**2)"),
                "theta_eq": 180 * u.deg,
            },
            independent_variables={"ri","rj","rk"},
        )

    @property
    def member_types(self):
        return self.__dict__.get("member_types_")

    @property
    def member_classes(self):
        return self.__dict__.get("member_classes_")


class VirtualPotentialType(ParametricPotential):
    """A description of non-bonded interactions between sites.

    This is a subclass of the gmso.core.Potential superclass.

    AtomType represents an atom type and includes the functional form
    describing its interactions and, optionally, other properties such as mass
    and charge.  This class inhereits from Potential, which stores the
    non-bonded interaction between atoms or sites. The functional form of the
    potential is stored as a `sympy` expression and the parameters, with units,
    are stored explicitly.
    """


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
        "", description="Description for the VirtualPotentialType", alias="description"
    )
    model_config = ConfigDict(
        alias_to_fields=dict(
            **ParametricPotential.model_config["alias_to_fields"],
            **{
        
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
        name="VirtualPotentialType",
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

        super(VirtualPotentialType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables,
            potential_expression=potential_expression,
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


    def clone(self, fast_copy=False):
        """Clone this AtomType, faster alternative to deepcopying."""
        return VirtualPotentialType(
            name=str(self.name),
            tags=self.tags,
            expression=None,
            parameters=None,
            independent_variables=None,
            potential_expression=self.potential_expression.clone(fast_copy),
            charge=u.unyt_quantity(self.charge.value, self.charge.units),
            atomclass=self.atomclass,
            doi=self.doi,
            overrides=(set(o for o in self.overrides) if self.overrides else None),
            description=self.description,
            definition=self.definition,
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
            and self.charge == other.charge
            and self.atomclass == other.atomclass
            and self.doi == other.doi
            and self.overrides == other.overrides
            and self.definition == other.definition
            and self.description == other.description
        )

    def _etree_attrib(self):
        attrib = super()._etree_attrib()
        if self.overrides == set():
            attrib.pop("overrides")
        charge = eval(attrib["charge"])
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



    @field_validator("charge_", mode="before")
    @classmethod
    def validate_charge(cls, charge):
        """Check to see that a charge is a unyt array of the right dimension."""
        if not isinstance(charge, u.unyt_array):
            warnings.warn(UNIT_WARNING_STRING.format("Charges", "elementary charge"))
            charge *= u.Unit("elementary_charge", registry=GMSO_UnitRegistry().reg)
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