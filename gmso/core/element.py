"""Representation of the chemical elements."""

import json
import warnings
from re import sub
from typing import Union

import numpy as np
import unyt as u
from pkg_resources import resource_filename
from pydantic import ConfigDict, Field, field_serializer

from gmso.abc.gmso_base import GMSOBase
from gmso.abc.serialization_utils import unyt_to_dict
from gmso.exceptions import GMSOError
from gmso.utils.misc import unyt_to_hashable

exported = [
    "element_by_mass",
    "element_by_symbol",
    "element_by_name",
    "element_by_smarts_string",
    "element_by_atomic_number",
    "element_by_atom_type",
    "Element",
]


class Element(GMSOBase):
    """Chemical element object

    Template to create a chemical element.
    Properties of the element instance are immutable.
    All known elements are pre-built and stored internally.
    """

    name: str = Field(..., description="Name of the element.")

    symbol: str = Field(..., description="Chemical symbol of the element.")

    atomic_number: int = Field(..., description="Atomic number of the element.")

    mass: u.unyt_quantity = Field(..., description="Mass of the element.")

    @field_serializer("mass")
    def serialize_mass(self, mass: Union[u.unyt_quantity, None]):
        if mass is None:
            return None
        else:
            return unyt_to_dict(mass)

    def __repr__(self):
        """Representation of the element."""
        return (
            f"<Element: {self.name}, symbol: {self.symbol}, "
            f"atomic number: {self.atomic_number}, mass: {self.mass.to('amu')}>"
        )

    def __eq__(self, other):
        if other is self:
            return True
        if not isinstance(other, Element):
            return False
        return (
            self.name == other.name
            and self.mass == other.mass
            and self.symbol == other.symbol
            and self.atomic_number == other.atomic_number
        )

    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)


def element_by_symbol(symbol, verbose=False):
    """Search for an element by its symbol.

    Look up an element from a list of known elements by symbol.
    Return None if no match found.

    Parameters
    ----------
    symbol : str
        Element symbol to look for, digits and spaces are removed before search.
    verbose : bool, optional, default=False
        If True, raise warnings if symbol has been trimmed before search.

    Returns
    -------
    matched_element : element.Element or None
        Return an element from the periodic table if the symbol is found,
        otherwise return None
    """
    symbol_trimmed = sub(r"[0-9 -]", "", symbol).capitalize()

    if symbol_trimmed != symbol and verbose:
        msg = (
            f"Numbers and spaces are not considered when searching by element symbol.\n"
            f"{symbol} became {symbol_trimmed}"
        )
        warnings.warn(msg)

    matched_element = symbol_dict.get(symbol_trimmed)
    return matched_element


def element_by_name(name, verbose=False):
    """Search for an element by its name.

    Look up an element from a list of known elements by name.
    Return None if no match found.

    Parameters
    ----------
    name : str
        Element name to look for, digits and spaces are removed before search.
    verbose : bool, optional, default=False
        If True, raise warnings if name has been trimmed before search.

    Returns
    -------
    matched_element : element.Element or None
        Return an element from the periodic table if the name is found,
        otherwise return None
    """
    name_trimmed = sub(r"[0-9 -]", "", name).lower()

    if name_trimmed != name and verbose:
        msg = (
            "Numbers and spaces are not considered when searching by element name. \n"
            f"{name} became {name_trimmed}"
        )
        warnings.warn(msg)

    matched_element = name_dict.get(name_trimmed)
    return matched_element


def element_by_atomic_number(atomic_number, verbose=False):
    """Search for an element by its atomic number.

    Look up an element from a list of known elements by atomic number.
    Return None if no match found.

    Parameters
    ----------
    atomic_number : int
        Element atomic number that need to look for
        if a string is provided, only numbers are considered during the search.
    verbose : bool, optional, default=False
        If True, raise warnings if atomic_number has been trimmed before search.

    Returns
    -------
    matched_element : element.Element
        Return an element from the periodic table if we find a match,
        otherwise raise GMSOError
    """
    if isinstance(atomic_number, str):
        atomic_number_trimmed = int(
            sub("[a-z -]", "", atomic_number.lower()).lstrip("0")
        )
        if str(atomic_number_trimmed) != atomic_number and verbose:
            msg = (
                f"Letters and spaces are not considered when searching by element atomic number. \n "
                f"{atomic_number} became {atomic_number_trimmed}"
            )
            warnings.warn(msg)
    else:
        atomic_number_trimmed = atomic_number
    matched_element = atomic_dict.get(atomic_number_trimmed)
    if matched_element is None:
        raise GMSOError(
            f"Failed to find an element with atomic number {atomic_number_trimmed}"
        )
    return matched_element


def element_by_mass(mass, exact=True, verbose=False):
    """Search for an element by its mass.

    Look up an element from a list of known elements by mass.
    If given mass is an int or a float, it will be convert to a unyt quantity (u.amu).
    Return None if no match found.

    Parameters
    ----------
    mass : int, float
        Element mass that need to look for, if a string is provided,
        only numbers are considered during the search.
        Mass unyt is assumed to be u.amu, unless specfied (which will be converted to u.amu).
    exact : bool, optional,  default=True
        This method can be used to search for an exact mass (up to the first decimal place)
        or search for an element  mass closest to the mass entered.
    verbose : bool, optional, default=False
        If True, raise warnings if mass has been trimmed before search.

    Returns
    -------
    matched_element : element.Element or None
        Return an element from the periodict table if we find a match,
        otherwise return None
    """
    if isinstance(mass, str):
        # Convert to float if a string is provided
        mass_trimmed = np.round(float(sub(r"[a-z -]", "", mass.lower())))
        if str(mass_trimmed) != mass and verbose:
            msg1 = (
                f"Letters and spaces are not considered when searching by element mass.\n"
                f"{mass} became {mass_trimmed}"
            )
            warnings.warn(msg1)
    elif isinstance(mass, u.unyt_quantity):
        # Convert to u.amu if a unyt_quantity is provided
        mass_trimmed = np.round(float(mass.to("amu")), 1)
    else:
        mass_trimmed = np.round(mass, 1)

    if exact:
        # Exact search mode
        matched_element = mass_dict.get(mass_trimmed)
    else:
        # Closest match mode
        mass_closest = min(
            mass_dict.keys(), key=lambda k: abs(k - mass_trimmed)
        )
        if verbose:
            msg2 = f"Closest mass to {mass_trimmed}: {mass_closest}"
            warnings.warn(msg2)
        matched_element = mass_dict.get(mass_closest)
    return matched_element


def element_by_smarts_string(smarts_string, verbose=False):
    """Search for an element by a given SMARTS string.

    Look up an element from a list of known elements by SMARTS string.
    Return None if no match found.

    Parameters
    ----------
    smarts_string : str
        SMARTS string representation of an atom type or its local chemical
        context. The Foyer SMARTS parser will be used to find the central atom
        and look up an Element. Note that this means some SMARTS grammar may
        not be parsed properly. For details, see
        https://github.com/mosdef-hub/foyer/issues/63
    verbose : bool, optional, default=False
        If True, raise warnings if smarts_string has been trimmed before search.

    Returns
    -------
    matched_element : element.Element
        Return an element from the periodic table if we find a match

    Raises
    ------
    GMSOError
        If no matching element is found for the provided smarts string
    """
    from gmso.utils.io import import_

    foyer = import_("foyer")
    SMARTS = foyer.smarts.SMARTS

    PARSER = SMARTS()

    symbols = PARSER.parse(smarts_string).iter_subtrees_topdown()

    first_symbol = None
    for symbol in symbols:
        if symbol.data == "atom_symbol":
            first_symbol = symbol.children[0]
            break

    matched_element = None
    if first_symbol is not None:
        matched_element = element_by_symbol(first_symbol)

    if matched_element is None:
        raise GMSOError(
            f"Failed to find an element from SMARTS string {smarts_string}. The "
            f"parser detected a central node with name {first_symbol}"
        )

    return matched_element


def element_by_atom_type(atom_type, verbose=False):
    """Search for an element by a given gmso AtomType object.

    Look up an element from a list of known elements by atom type.
    Return None if no match is found.

    Parameters
    ----------
    atom_type : gmso.core.atom_type.AtomType
        AtomType object to be parsed for element information. Attributes are
        looked up in the order of mass, name, and finally definition (the
        SMARTS string).  Because of the loose structure of this class, a
        successful lookup is not guaranteed.
    verbose : bool, optional, default=False
        If True, raise warnings if atom_type has been trimmed before search.

    Returns
    -------
    matched_element : element.Element or None
        Return an element from the periodict table if we find a match,
        otherwise return None

    """
    matched_element = None

    if matched_element is None and atom_type.mass:
        matched_element = element_by_mass(
            atom_type.mass, exact=False, verbose=verbose
        )
    if matched_element is None and atom_type.name:
        matched_element = element_by_symbol(atom_type.name, verbose=verbose)
    if matched_element is None and atom_type.definition:
        matched_element = element_by_smarts_string(
            atom_type.definition, verbose=verbose
        )

    if matched_element is None:
        raise GMSOError(
            f"Failed to find an element from atom type"
            "{atom_type} with "
            "properties mass: {atom_type.mass}, name:"
            "{atom_type.name}, and "
            "definition: {atom_type.definition}"
        )

    return matched_element


# Get the JSON file from ele package for a standard representation
elements_json_loc = resource_filename("ele", "lib/elements.json")
elements_dict = None
elements = []
with open(elements_json_loc, "r") as el_json_file:
    elements_dict = json.load(el_json_file)
    elements_dict = {int(key): value for key, value in elements_dict.items()}


for atom_number, element_properties in elements_dict.items():
    assert atom_number == element_properties["atomic number"]
    elements.append(
        Element(
            atomic_number=element_properties["atomic number"],
            name=element_properties["name"],
            symbol=element_properties["symbol"],
            mass=element_properties["mass"] * u.amu,
        )
    )

elements.sort(key=lambda el: el.atomic_number)

symbol_dict = {element.symbol: element for element in elements}
name_dict = {element.name: element for element in elements}
atomic_dict = {element.atomic_number: element for element in elements}
mass_dict = {
    np.round(float(element.mass.to("amu")), 1): element for element in elements
}
element_by_capitalized_names = {
    element.name.capitalize(): element for element in elements
}


def __getattr__(name):
    """Access an element by the object attribute."""
    if name in exported:
        return globals()[name]
    elif name in element_by_capitalized_names:
        return element_by_capitalized_names[name]
    else:
        raise AttributeError(f"Module {__name__} has no attribute {name}")
