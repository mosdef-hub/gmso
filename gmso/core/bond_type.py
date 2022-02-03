"""Defined interactions between two bonded partners."""
from typing import Optional, Tuple

import unyt as u
from pydantic import Field

from gmso.core.parametric_potential import ParametricPotential
from gmso.utils._constants import BOND_TYPE_DICT
from gmso.utils.expression import PotentialExpression


class BondType(ParametricPotential):
    __base_doc__ = """A descripton of the interaction between 2 bonded partners.

    This is a subclass of the gmso.core.Potential superclass.

    BondType represents a bond type and includes the functional form describing
    its interactions. The functional form of the potential is stored as a
    `sympy` expression and the parameters, with units, are stored explicitly.
    The AtomTypes that are used to define the bond type are stored as
    `member_types`.

    Notes
    ----
    Inherits many functions from gmso.ParametricPotential:
        __eq__, _validate functions
    """

    member_types_: Optional[Tuple[str, str]] = Field(
        None,
        description="List-like of of gmso.AtomType.name "
        "defining the members of this bond type",
    )

    member_classes_: Optional[Tuple[str, str]] = Field(
        None,
        description="List-like of of gmso.AtomType.atomclass "
        "defining the members of this bond type",
    )

    def __init__(
        self,
        name="BondType",
        expression=None,
        parameters=None,
        independent_variables=None,
        potential_expression=None,
        member_types=None,
        member_classes=None,
        topology=None,
        tags=None,
    ):
        super(BondType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables,
            potential_expression=potential_expression,
            topology=topology,
            member_types=member_types,
            member_classes=member_classes,
            set_ref=BOND_TYPE_DICT,
            tags=tags,
        )

    @property
    def member_types(self):
        """Return the members involved in this bondtype."""
        return self.__dict__.get("member_types_")

    @property
    def member_classes(self):
        return self.__dict__.get("member_classes_")

    @staticmethod
    def _default_potential_expr():
        return PotentialExpression(
            expression="0.5 * k * (r-r_eq)**2",
            independent_variables={"r"},
            parameters={
                "k": 1000 * u.Unit("kJ / (nm**2)"),
                "r_eq": 0.14 * u.nm,
            },
        )

    class Config:
        """Pydantic configuration for class attributes."""

        fields = {
            "member_types_": "member_types",
            "member_classes_": "member_classes",
        }

        alias_to_fields = {
            "member_types": "member_types_",
            "member_classes": "member_classes_",
        }
