"""Abstract representation of a Potential object."""

from abc import abstractmethod
from typing import Any, Dict, Iterator, List

import unyt as u
from pydantic import ConfigDict, Field, field_serializer, field_validator

from gmso.abc.gmso_base import GMSOBase
from gmso.abc.serialization_utils import unyt_to_dict
from gmso.utils.expression import PotentialExpression


class AbstractPotential(GMSOBase):
    __base_doc__ = """An abstract potential class.

    AbstractPotential stores a general interaction between components of a chemical
    topology that can be specified by a mathematical expression. The functional
    form of the potential is stored as a `sympy` expression and the parameters
    are stored explicitly. This class is agnostic to the instantiation of the
    potential, which can be e.g. a non-bonded potential, a bonded potential, an
    angle potential, a dihedral potential, etc. and is designed to be inherited
    by classes that represent these potentials.
    """

    name_: str = Field(
        "",
        description="The name of the potential. Defaults to class name",
        alias="name",
    )

    potential_expression_: PotentialExpression = Field(
        PotentialExpression(expression="a*x+b", independent_variables={"x"}),
        description="The mathematical expression for the potential",
        alias="potential_expression",
    )

    tags_: Dict[str, Any] = Field(
        {},
        description="Tags associated with the potential",
        alias="tags",
    )
    model_config = ConfigDict(
        alias_to_fields={
            "name": "name_",
            "potential_expression": "potential_expression_",
            "tags": "tags_",
        }
    )

    def __init__(
        self,
        name="Potential",
        expression="a*x+b",
        independent_variables=None,
        potential_expression=None,
        **kwargs,
    ):
        if potential_expression is None:
            if expression is None:
                expression = "a*x+b"

            if independent_variables is None:
                independent_variables = {"x"}

            potential_expression = PotentialExpression(
                expression=expression,
                independent_variables=independent_variables,
                parameters=None,
            )

        if not kwargs.get("tags"):
            kwargs["tags"] = {}

        super().__init__(
            name=name, potential_expression=potential_expression, **kwargs
        )

    @property
    def name(self):
        """The name of the potential."""
        return self.__dict__.get("name_")

    @property
    def independent_variables(self):
        """Optional[Union[set, str]]\n\tThe independent variables in the `Potential`'s expression."""
        return self.potential_expression.independent_variables

    @property
    def expression(self):
        """Optional[Union[str, sympy.Expr]]\n\tThe mathematical expression of the functional form of the potential."""
        return self.potential_expression.expression

    @property
    def potential_expression(self):
        """Return the functional form of the potential."""
        return self.__dict__.get("potential_expression_")

    @property
    def tags(self):
        return self.__dict__.get("tags_")

    @property
    def tag_names(self) -> List[str]:
        return list(self.__dict__.get("tags_"))

    @property
    def tag_names_iter(self) -> Iterator[str]:
        return iter(self.__dict__.get("tags_"))

    @field_serializer("potential_expression_")
    def serialize_expression(self, potential_expression_: PotentialExpression):
        expr = str(potential_expression_.expression)
        ind = sorted(
            list(
                str(ind) for ind in potential_expression_.independent_variables
            )
        )
        params = {
            param: unyt_to_dict(val)
            for param, val in potential_expression_.parameters.items()
        }
        return {
            "expression": expr,
            "independent_variables": ind,
            "parameters": params,
        }

    @field_serializer("tags_")
    def serialize_tags(self, tags_):
        return_dict = dict()
        for key, val in tags_.items():
            if isinstance(val, u.unyt_array):
                return_dict[key] = unyt_to_dict(val)
            else:
                return_dict[key] = val
        return return_dict

    def add_tag(self, tag: str, value: Any, overwrite=True) -> None:
        """Add metadata for a particular tag"""
        if self.tags.get(tag) and not overwrite:
            raise ValueError(
                f"Tag {tag} already exists. "
                f"Please use overwrite=True to overwrite"
            )
        self.tags[tag] = value

    def get_tag(self, tag: str, throw=False) -> Any:
        """Get value of a particular tag"""
        if throw:
            return self.tags[tag]
        else:
            return self.tags.get(tag)

    def delete_tag(self, tag: str) -> None:
        del self.tags[tag]

    def pop_tag(self, tag: str) -> Any:
        return self.tags.pop(tag, None)

    @field_validator("potential_expression_", mode="before")
    @classmethod
    def validate_potential_expression(cls, v):
        if isinstance(v, dict):
            v = PotentialExpression(**v)
        return v

    @abstractmethod
    def set_expression(self):
        """Set the functional form of the expression."""
        raise NotImplementedError

    def __setattr__(self, key: Any, value: Any) -> None:
        """Set attributes of the potential."""
        if key == "expression":
            self.potential_expression.expression = value
        elif key == "independent_variables":
            self.potential_expression.independent_variables = value
        elif key == "set_ref_":
            return
        else:
            super().__setattr__(key, value)

    def __repr__(self):
        """Return a formatted representation of the potential."""
        desc = (
            f"<{self.__class__.__name__} {self.name},\n "
            f"expression: {self.expression},\n "
            f"id: {id(self)}>"
        )
        return desc

    def __str__(self):
        """Return a string representation of the potential."""
        return (
            f"<{self.__class__.__name__} {self.name}, "
            f"expression: {self.expression}, "
            f"id: {id(self)}>"
        )
