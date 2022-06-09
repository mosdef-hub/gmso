from typing import Optional, Tuple

import unyt as u
from pydantic import Field

from gmso.core.parametric_potential import ParametricPotential
from gmso.utils._constants import ANGLE_TYPE_DICT
from gmso.utils.expression import PotentialExpression


class AngleType(ParametricPotential):
    __base_doc__ = """A descripton of the interaction between 3 bonded partners.

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
    )

    member_classes_: Optional[Tuple[str, str, str]] = Field(
        None,
        description="List-like of gmso.AtomType.atomclass "
        "defining the members of this angle type",
    )

    def __init__(
        self,
        name="AngleType",
        expression=None,
        parameters=None,
        independent_variables=None,
        potential_expression=None,
        member_types=None,
        member_classes=None,
        topology=None,
        tags=None,
    ):

        super(AngleType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables,
            potential_expression=potential_expression,
            topology=topology,
            member_types=member_types,
            member_classes=member_classes,
            set_ref=ANGLE_TYPE_DICT,
            tags=tags,
        )

    @staticmethod
    def _default_potential_expr():
        return PotentialExpression(
            expression="0.5 * k * (theta-theta_eq)**2",
            parameters={
                "k": 1000 * u.Unit("kJ / (deg**2)"),
                "theta_eq": 180 * u.deg,
            },
            independent_variables={"theta"},
        )

    @property
    def member_types(self):
        return self.__dict__.get("member_types_")

    @property
    def member_classes(self):
        return self.__dict__.get("member_classes_")

    class Config:
        fields = {
            "member_types_": "member_types",
            "member_classes_": "member_classes",
        }

        alias_to_fields = {
            "member_types": "member_types_",
            "member_classes": "member_classes_",
        }
