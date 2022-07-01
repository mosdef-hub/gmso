from copy import copy, deepcopy
from typing import Any, Union

import unyt as u

from gmso.abc.abstract_potential import AbstractPotential
from gmso.utils.expression import PotentialExpression


class ParametricPotential(AbstractPotential):
    __base_doc__ = """A parametric potential class.

    Potential stores a general interaction between components of a chemical
    topology that can be specified by a mathematical expression. The functional
    form of the potential is stored as a `sympy` expression and the parameters
    are stored explicitly. This class is agnostic to the instantiation of the
    potential, which can be e.g. a non-bonded potential, a bonded potential, an
    angle potential, a dihedral potential, etc. and is designed to be inherited
    by classes that represent these potentials.
    """

    def __init__(
        self,
        name="ParametricPotential",
        expression=None,
        parameters=None,
        potential_expression=None,
        independent_variables=None,
        **kwargs,
    ):
        if potential_expression is not None and (
            expression is not None
            or independent_variables is not None
            or parameters is not None
        ):
            raise ValueError(
                "When using potential expressions "
                "please do not provide arguments for "
                "expression, independent_variables or parameters."
            )

        if potential_expression is None:
            _potential_expression = self._get_expression(
                expression, parameters, independent_variables
            )
        else:
            _potential_expression = potential_expression

        super().__init__(
            name=name,
            potential_expression=_potential_expression,
            **kwargs,
        )

    def _get_expression(self, expression, parameters, indep_vars):
        args = (expression, parameters, indep_vars)
        all_provided = tuple(1 if param is not None else 0 for param in args)

        if sum(all_provided) == 0:
            return self._default_potential_expr()
        elif sum(all_provided) < 3:
            raise ValueError(
                "When using keyword arguments `expression`, "
                "`independent_variables` and `parameters` for "
                "a potential, you are expected to provide all the values "
                "or none of them to use defaults. However, you provided "
                "the following and there's not enough information to form "
                "a set of expression, idependent_variables and parameters.\n"
                f"expression: {expression}\n"
                f"parameters: {parameters}\n"
                f"independent_variables: {indep_vars}\n"
            )
        else:
            return PotentialExpression(
                expression=expression,
                independent_variables=indep_vars,
                parameters=parameters,
            )

    @staticmethod
    def _default_potential_expr():
        return PotentialExpression(
            expression="a*x+b",
            parameters={"a": 1.0 * u.dimensionless, "b": 1.0 * u.dimensionless},
            independent_variables={"x"},
        )

    @property
    def parameters(self):
        """Optional[dict]\n\tThe parameters of the `Potential` expression and their corresponding values, as `unyt` quantities"""
        return self.potential_expression_.parameters

    def __setattr__(self, key: Any, value: Any) -> None:
        """Set the attributes of the potential."""
        if key == "parameters":
            self.potential_expression_.parameters = value
        else:
            super().__setattr__(key, value)

    def set_expression(
        self, expression=None, parameters=None, independent_variables=None
    ):
        """Set the expression, parameters, and independent variables for this potential.

        Parameters
        ----------
        expression: sympy.Expression or string
            The mathematical expression corresponding to the potential
            If None, the expression remains unchanged
        parameters: dict
            {parameter: value} in the expression
            If None, the parameters remain unchanged

        Notes
        -----
        Be aware of the symbols used in the `expression` and `parameters`.
        If unnecessary parameters are supplied, an error is thrown.
        If only a subset of the parameters are supplied, they are updated
        while the non-passed parameters default to the existing values
        """
        self.potential_expression_.set(
            expression=expression,
            independent_variables=independent_variables,
            parameters=parameters,
        )

    def dict(
        self,
        *,
        include: Union["AbstractSetIntStr", "MappingIntStrAny"] = None,
        exclude: Union["AbstractSetIntStr", "MappingIntStrAny"] = None,
        by_alias: bool = False,
        skip_defaults: bool = None,
        exclude_unset: bool = False,
        exclude_defaults: bool = False,
        exclude_none: bool = False,
    ) -> dict:
        if exclude is None:
            exclude = set()
        if isinstance(exclude, dict):
            exclude = set(exclude)

        exclude = exclude.union({"topology_", "set_ref_"})

        return super().dict(
            include=include,
            exclude=exclude,
            by_alias=True,
            skip_defaults=skip_defaults,
            exclude_unset=exclude_unset,
            exclude_defaults=exclude_defaults,
            exclude_none=exclude_none,
        )

    def __eq__(self, other):
        if other is self:
            return True
        if not isinstance(other, type(self)):
            return False
        return (
            self.expression == other.expression
            and self.independent_variables == other.independent_variables
            and self.name == other.name
            and self.parameters == other.parameters
        )

    def get_parameters(self, copy=False):
        """Return parameters for this ParametricPotential."""
        if copy:
            params = {
                k: u.unyt_quantity(v.value, v.units)
                for k, v in self.parameters.items()
            }
        else:
            params = self.parameters

        return params

    def clone(self, fast_copy=False):
        """Clone this parametric potential, faster alternative to deepcopying."""
        Creator = self.__class__
        kwargs = {"tags": deepcopy(self.tags_)}
        if hasattr(self, "member_classes"):
            kwargs["member_classes"] = (
                copy(self.member_classes) if self.member_classes else None
            )

        if hasattr(self, "member_types"):
            kwargs["member_types"] = (
                copy(self.member_types) if self.member_types else None
            )

        return Creator(
            name=self.name,
            potential_expression=self.potential_expression_.clone(fast_copy),
            **kwargs,
        )

    @classmethod
    def from_template(cls, potential_template, parameters, name=None, **kwargs):
        """Create a potential object from the potential_template.

        Parameters
        ----------
        potential_template : gmso.lib.potential_templates.PotentialTemplate,
            The potential template object
        parameters : dict,
            The parameters of the potential object to create
        name: str, default=None,
            The new name for the created parametric potential, defaults to the
            template name.
        **kwargs: dict
            The remaining keyword arguments to the Parametric potential's constructor.

        Returns
        -------
        gmso.ParametricPotential
            The potential object created

        Raises
        ------
        GMSOError
            If potential_template is not of instance PotentialTemplate
        """
        from gmso.lib.potential_templates import PotentialTemplate

        if not isinstance(potential_template, PotentialTemplate):
            raise TypeError(
                f"Object {potential_template} of type {type(potential_template)} is not an instance of "
                f"PotentialTemplate."
            )

        potential_template.assert_can_parameterize_with(parameters)
        new_expression = PotentialExpression.from_non_parametric(
            potential_template.potential_expression, parameters, valid=True
        )
        return cls(
            name=name or potential_template.name,
            potential_expression=new_expression,
            **kwargs,
        )

    def __repr__(self):
        """Return formatted representation of the potential."""
        desc = super().__repr__()
        member_types = (
            lambda x: x.member_types if hasattr(x, "member_types") else ""
        )
        desc = desc.replace(
            ">",
            f", \n parameters: {self.parameters},\n"
            f"member types: {member_types(self)}>",
        )
        return desc

    class Config:
        """Pydantic configuration class."""

        fields = {"topology_": "topology", "set_ref_": "set_ref"}
        alias_to_fields = {"topology": "topology_"}
        validate_assignment = True
