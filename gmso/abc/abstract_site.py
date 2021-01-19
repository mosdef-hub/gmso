import warnings
from typing import Union, Sequence, Any, TypeVar, ClassVar

import numpy as np
import unyt as u
from unyt.exceptions import InvalidUnitOperation

from pydantic import validator, Field

from gmso.abc.gmso_base import GMSOBase
from gmso.exceptions import GMSOError

PositionType = Union[Sequence[float], np.ndarray, u.unyt_array]
SiteT = TypeVar('SiteT', bound='Site')

BASE_DOC_ATTR = '__base_doc__'
FIELDS_IN_DOCSTRING = 'alias_to_fields'


def default_position():
    return u.unyt_array([np.nan] * 3, u.nm)


class Site(GMSOBase):
    __base_doc__: ClassVar[str] = """An interaction site object in the topology hierarchy.

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
        '',
        description='Name of the site, defaults to class name',
    )

    label_: str = Field(
        '',
        description='Label to be assigned to the site'
    )

    position_: PositionType = Field(
        default_factory=default_position,
        description='The 3D Cartesian coordinates of the position of the site'
    )

    @property
    def name(self) -> str:
        return self.__dict__.get('name_')

    @property
    def position(self) -> u.unyt_array:
        return self.__dict__.get('position_')

    @property
    def label(self) -> str:
        return self.__dict__.get('label_')

    def __repr__(self):
        return f'<{self.__class__.__name__}, id {id(self)}>'

    @validator('position_')
    def is_valid_position(cls, position):
        """Validator for attribute position"""
        if position is None:
            return u.unyt_array([np.nan] * 3, u.nm)

        if not isinstance(position, u.unyt_array):
            try:
                position *= u.nm
            except InvalidUnitOperation as e:
                raise GMSOError(f'Converting object of type {type(position)} failed with following error: {e}')
            warnings.warn('Positions are assumed to be in nm')

        try:
            position = np.reshape(position, newshape=(3,), order='C')
            position.convert_to_units(u.nm)
        except ValueError:
            raise ValueError(f'Position of shape {position.shape} is not valid. '
                             'Accepted values: (a.) list-like of length 3'
                             '(b.) np.array or unyt.unyt_array of shape (3,)')

        return position

    @validator('name_', pre=True, always=True)
    def inject_name(cls, value):
        if value == '' or value is None:
            return cls.__name__
        else:
            return value

    @classmethod
    def __new__(cls, *args: Any, **kwargs: Any) -> SiteT:
        if cls is Site:
            raise TypeError('Cannot instantiate abstract class of type Site')
        else:
            return object.__new__(cls)

    class Config:
        arbitrary_types_allowed = True

        fields = {
            'name_': 'name',
            'position_': 'position',
            'label_': 'label'
        }

        alias_to_fields = {
            'name': 'name_',
            'position': 'position_',
            'label': 'label_'
        }

        validate_assignment = True
