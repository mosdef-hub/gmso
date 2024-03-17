from typing import Optional, Tuple

import unyt as u
from pydantic import ConfigDict, Field

from gmso.core.parametric_potential import ParametricPotential
from gmso.utils.expression import PotentialExpression


class PairPotentialType(ParametricPotential):
    """A description of custom pairwise potential between 2 AtomTypes that does not follow combination rule.

    This is a subclass of the gmso.core.ParametricPotential superclass.

    PairPotentialType represents a type of pairwise potential between two
    Atomtypes that does not follow a specific combination rule, and includes the functional
    form describing its interactions. The functional form of the potential is
    stored as a `sympy` expression and the parameters, with units, are stored
    explicitly.  The AtomTypes that are used to define the dihedral type are
    stored as `member_types`.


    Notes
    -----
    Inherits many functions from gmso.ParametricPotential:
        __eq__, _validate functions
    """

    member_types_: Optional[Tuple[str, str]] = Field(
        None,
        description="List-like of strs, referring to gmso.Atomtype.name or gmso.Atomtype.atomclass, "
        "defining the members of this pair potential type",
        alias="member_types",
    )
    model_config = ConfigDict(
        alias_to_fields=dict(
            **ParametricPotential.model_config["alias_to_fields"],
            **{
                "member_types": "member_types_",
            },
        ),
    )

    def __init__(
        self,
        name="PairPotentialType",
        expression=None,
        parameters=None,
        independent_variables=None,
        potential_expression=None,
        member_types=None,
        tags=None,
    ):
        super(PairPotentialType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables,
            member_types=member_types,
            potential_expression=potential_expression,
            tags=tags,
        )

    @property
    def member_types(self):
        return self.__dict__.get("member_types_")

    @staticmethod
    def _default_potential_expr():
        return PotentialExpression(
            expression="4 * eps * ((sigma / r)**12 - (sigma / r)**6)",
            independent_variables={"r"},
            parameters={"eps": 1 * u.Unit("kJ / mol"), "sigma": 1 * u.nm},
        )
