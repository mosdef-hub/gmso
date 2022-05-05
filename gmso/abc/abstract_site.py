"""Basic interaction site in GMSO that all other sites will derive from."""
import warnings
from typing import Any, ClassVar, Optional, Sequence, TypeVar, Union

import numpy as np
import unyt as u
from pydantic import Field, StrictInt, StrictStr, validator
from unyt.exceptions import InvalidUnitOperation

from gmso.abc.gmso_base import GMSOBase
from gmso.exceptions import GMSOError

PositionType = Union[Sequence[float], np.ndarray, u.unyt_array]
SiteT = TypeVar("SiteT", bound="Site")

BASE_DOC_ATTR = "__base_doc__"
FIELDS_IN_DOCSTRING = "alias_to_fields"


def default_position():
    return u.unyt_array([np.nan] * 3, u.nm)


class Site(GMSOBase):
    __iterable_attributes__: ClassVar[set] = {
        "label",
        "group",
        "molecule",
        "residue_name",
        "residue_number",
    }

    __base_doc__: ClassVar[
        str
    ] = """An interaction site object in the topology hierarchy.

    Site is the object that represents any general interaction site in a molecular simulation.
    Sites have been designed to be as general as possible, making no assumptions about representing atoms or beads, or
    having mass or charge. That is, a Site can represent an atom in an atomistic system,
    a bead in a coarse-grained system, and much more.

    Notes
    -----
    The label attribute for a site takes its meaning when used with some sort of container (like topology)
    such that a label for a site can then be used to group sites together. The rules for defining a site label
    and their meaning the responsibility of the container where the sites will reside.
    """

    name_: str = Field(
        "",
        description="Name of the site, defaults to class name",
    )

    label_: str = Field("", description="Label to be assigned to the site")

    group_: str = Field(
        "DefaultGroup", description="Molecule Group label for the site"
    )

    molecule_name_: Optional[StrictStr] = Field(
        "DefaultMolecule", description="Molecule label for the site"
    )

    molecule_number_: Optional[StrictInt] = Field(
        None, description="Molecule number for the site"
    )

    residue_name_: Optional[StrictStr] = Field(
        None, description="Residue label for the site"
    )

    residue_number_: Optional[StrictInt] = Field(
        None, description="Residue number for the site"
    )

    position_: PositionType = Field(
        default_factory=default_position,
        description="The 3D Cartesian coordinates of the position of the site",
    )

    @property
    def name(self) -> str:
        """Return the name of the site."""
        return self.__dict__.get("name_")

    @property
    def position(self) -> u.unyt_array:
        """Return the 3D Cartesian coordinates of the site."""
        return self.__dict__.get("position_")

    @property
    def label(self) -> str:
        """Return the label assigned to the site."""
        return self.__dict__.get("label_")

    @property
    def group(self) -> str:
        """Return the group of the site."""
        return self.__dict__.get("group_")

    @property
    def molecule_name(self) -> str:
        """Return the molecule name of the site."""
        return self.__dict__.get("molecule_name_")

    @property
    def molecule_number(self) -> str:
        """Return the molecule number of the site."""
        return self.__dict__.get("molecule_number_")

    @property
    def residue_name(self):
        """Return the residue name assigned to the site."""
        return self.__dict__.get("residue_name_")

    @property
    def residue_number(self):
        """Return the reside number assigned to the site."""
        return self.__dict__.get("residue_number_")

    def __repr__(self):
        """Return the formatted representation of the site."""
        return (
            f"<{self.__class__.__name__} {self.name},\n "
            f"position: {self.position},\n "
            f"label: {self.label if self.label else None},\n "
            f"id: {id(self)}>"
        )

    def __str__(self):
        """Return the string representation of the site."""
        return (
            f"<{self.__class__.__name__} {self.name}, "
            f"label: {self.label if self.label else None} id: {id(self)}>"
        )

    @validator("position_")
    def is_valid_position(cls, position):
        """Validate attribute position."""
        if position is None:
            return u.unyt_array([np.nan] * 3, u.nm)

        if not isinstance(position, u.unyt_array):
            try:
                position *= u.nm
            except InvalidUnitOperation as e:
                raise GMSOError(
                    f"Converting object of type {type(position)} failed with following error: {e}"
                )
            warnings.warn("Positions are assumed to be in nm")

        try:
            position = np.reshape(position, newshape=(3,), order="C")
            position.convert_to_units(u.nm)
        except ValueError:
            raise ValueError(
                f"Position of shape {position.shape} is not valid. "
                "Accepted values: (a.) list-like of length 3"
                "(b.) np.array or unyt.unyt_array of shape (3,)"
            )

        return position

    @validator("name_", pre=True, always=True)
    def inject_name(cls, value):
        if value == "" or value is None:
            return cls.__name__
        else:
            return value

    @classmethod
    def __new__(cls, *args: Any, **kwargs: Any) -> SiteT:
        if cls is Site:
            raise TypeError("Cannot instantiate abstract class of type Site")
        else:
            return object.__new__(cls)

    class Config:
        """Pydantic configuration for site objects."""

        arbitrary_types_allowed = True

        fields = {
            "name_": "name",
            "position_": "position",
            "label_": "label",
            "group_": "group",
            "molecule_name_": "molecule_name",
            "molecule_number_": "molecule_number",
            "residue_name_": "residue_name",
            "residue_number_": "residue_number",
        }

        alias_to_fields = {
            "name": "name_",
            "position": "position_",
            "label": "label_",
            "group": "group_",
            "molecule_name": "molecule_name_",
            "molecule_number": "molecule_number_",
            "residue_name": "residue_name_",
            "residue_number": "residue_number_",
        }

        validate_assignment = True
