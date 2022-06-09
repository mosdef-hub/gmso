"""Module supporting template potential objects."""
import json
from pathlib import Path

from gmso.abc.abstract_potential import AbstractPotential
from gmso.exceptions import GMSOError
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

    def __init__(
        self,
        name="PotentialTemplate",
        expression="4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
        independent_variables="r",
        potential_expression=None,
        template=True,
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
            name=name, potential_expression=_potential_expression
        )

    def set_expression(self, *args, **kwargs):
        """Set the expression of the PotentialTemplate."""
        raise NotImplementedError

    class Config:
        """Pydantic configuration for potential template."""

        allow_mutation = False


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
