"""Module supporting template potential objects."""

import json
from pathlib import Path
from typing import Dict

import sympy
import unyt as u
from pydantic import ConfigDict, Field, field_validator

from gmso.abc.abstract_potential import AbstractPotential
from gmso.exceptions import (
    GMSOError,
    MissingParameterError,
    UnknownParameterError,
)
from gmso.utils.expression import PotentialExpression
from gmso.utils.singleton import Singleton

POTENTIAL_JSONS = list(Path(__file__).parent.glob("jsons/*.json"))
JSON_DIR = Path.joinpath(Path(__file__).parent, "jsons")


def _verify_potential_template_keys(_dict, name):
    """Verify the potential template is properly formatted."""
    assert (
        "name" in _dict
    ), f"Key name not found in the potential template {name}.json"
    assert (
        "expression" in _dict
    ), f"Key expression not found in the potential template {name}.json"
    assert (
        "independent_variables" in _dict
    ), f"Key independent_variables not found in the potential template {name}.json"
    if str(name) != _dict["name"]:
        raise GMSOError(
            f'Mismatch between Potential name {name} and {_dict["name"]}'
        )


def _load_template_json(item, json_dir=JSON_DIR):
    """Return dictionary representation of PotentialTemplate from JSON."""
    with json_dir.joinpath(f"{item}.json").open("r") as json_file:
        potential_dict = json.load(json_file)
        _verify_potential_template_keys(potential_dict, item)
        return potential_dict


class PotentialTemplate(AbstractPotential):
    """Template for potential objects to be re-used."""

    expected_parameters_dimensions_: Dict[str, sympy.Expr] = Field(
        ...,
        description="The expected dimensions for parameters.",
        alias="expected_parameters_dimensions",
    )

    model_config = ConfigDict(
        frozen=True,
        alias_to_fields=dict(
            **AbstractPotential.model_config["alias_to_fields"],
            **{
                "expected_parameters_dimensions": "expected_parameters_dimensions_",
            },
        ),
    )

    def __init__(
        self,
        name="PotentialTemplate",
        expression="4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
        independent_variables="r",
        potential_expression=None,
        expected_parameters_dimensions=None,
    ):
        if not isinstance(independent_variables, set):
            independent_variables = set(independent_variables.split(","))

        if potential_expression is None:
            _potential_expression = PotentialExpression(
                expression=expression,
                independent_variables=independent_variables,
            )
        else:
            _potential_expression = potential_expression

        super(PotentialTemplate, self).__init__(
            name=name,
            potential_expression=_potential_expression,
            expected_parameters_dimensions=expected_parameters_dimensions,
        )

    @field_validator("expected_parameters_dimensions_", mode="before")
    def validate_expected_parameters(cls, dim_dict):
        """Validate the expected parameters and dimensions for this template."""
        if not isinstance(dim_dict, Dict):
            raise TypeError(
                f"Expected expected_parameters_dimensions to be a "
                f"dictionary but found {type(dim_dict)}"
            )
        for param_name, dim in dim_dict.items():
            if not isinstance(dim, sympy.Expr):
                try:
                    dimension = getattr(u.dimensions, dim)
                except AttributeError:
                    dimension_expr = sympy.sympify(dim)
                    subs = (
                        (symbol, getattr(u.dimensions, str(symbol)))
                        for symbol in dimension_expr.free_symbols
                    )
                    dimension = dimension_expr.subs(subs)
                dim_dict[param_name] = dimension
        return dim_dict

    @property
    def expected_parameters_dimensions(self):
        """Return the expected dimensions of the parameters for this template."""
        return self.__dict__.get("expected_parameters_dimensions_")

    def set_expression(self, *args, **kwargs):
        """Set the expression of the PotentialTemplate."""
        raise NotImplementedError

    def assert_can_parameterize_with(
        self, parameters: Dict[str, u.unyt_quantity]
    ) -> None:
        """Assert that a ParametricPotential can be instantiated from this template and provided parameters."""
        if not isinstance(parameters, dict):
            raise TypeError("Provided `parameters` is not a dictionary.")

        for param_name in self.expected_parameters_dimensions:
            if param_name not in parameters:
                raise MissingParameterError(param_name, list(parameters))

        for param_name, param_value in parameters.items():
            quantity = param_value
            if not (isinstance(param_value, u.unyt_array)):
                raise ValueError(f"Parameter {param_name} lacks a unit.")

            if param_name not in self.expected_parameters_dimensions:
                raise UnknownParameterError(
                    param_name, list(self.expected_parameters_dimensions)
                )
            expected_param_dimension = self.expected_parameters_dimensions[
                param_name
            ]
            param_dimension = quantity.units.dimensions
            if param_dimension != expected_param_dimension:
                if expected_param_dimension == 1:
                    expected_param_dimension = "dimensionless"
                if param_dimension == 1:
                    param_dimension = "dimensionless"

                raise AssertionError(
                    f"Expected parameter {param_name} to have "
                    f"dimension {expected_param_dimension} but found {param_dimension}. "
                    f"So, a {self.__class__.__name__} cannot be instantiated using the provided "
                    f"parameters: {parameters}"
                )


class PotentialTemplateLibrary(Singleton):
    """A singleton collection of all the potential templates."""

    def __init__(self):
        try:
            self.json_refs
        except AttributeError:
            self.json_refs = POTENTIAL_JSONS
            potential_names = [pot_json.name for pot_json in POTENTIAL_JSONS]
            self._ref_dict = {
                potential.replace(".json", ""): None
                for potential in potential_names
            }
            self._user_ref_dict = {}

    def get_available_template_names(self):
        """Return all available potential templates available to gmso."""
        return tuple(self._ref_dict.keys())

    def __getitem__(self, item):
        """Return template if it exists."""
        if item not in self._ref_dict:
            raise KeyError(f"Potential Template {item} not found.")

        return self._load_from_json(item)

    def _load_from_json(self, item):
        """Return template from JSON."""
        if self._ref_dict[item] is None:
            potential_dict = _load_template_json(item, JSON_DIR)
            self._ref_dict[item] = PotentialTemplate(**potential_dict)
        return self._ref_dict[item]
