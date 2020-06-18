import warnings
from abc import ABC
from typing import Union, Sequence, Optional, Any, TypeVar, ClassVar

import numpy as np
import unyt as u

from pydantic import BaseModel, validator, root_validator, Field

from gmso.abc.auto_doc import apply_docs


PositionType = Union[Sequence[float], np.ndarray]
SiteT = TypeVar('SiteT', bound='Site')

BASE_DOC_ATTR = '__base_doc__'
FIELDS_IN_DOCSTRING = 'alias_to_fields'


class Site(BaseModel, ABC):

    __base_doc__: ClassVar[str] = """An interaction site object in the topology hierarchy.
    
    Site is the object that represents any general interaction site in a molecular simulation.
    Sites have been designed to be as general as possible, making no assumptions about representing atoms or beads, or
    having mass or charge. That is, a Site can represent an atom in an atomistic system,
    a bead in a coarse-grained system, and much more.
    
    Notes
    -----
    The label attribute for a site takes its meaning when used with some sort of container (like topology)
    such that a label for a site can then be used to group sites together. The rules for defining a site label
    and their meaning is left upto the container where the sites will reside
    """

    __docs_generated__: ClassVar[bool] = False

    name_: str = Field(
        '',
        description='Name of the site, defaults to class name',
    )

    label_: str = Field(
        '',
        description='Label to be assigned for the site'
    )

    position_: Optional[PositionType] = Field(
        None,
        description='The Cartesian coordinates of the position of the site'
    )

    @property
    def name(self) -> str:
        return self.__dict__.get('name_')

    @property
    def position(self) -> PositionType:
        return self.__dict__.get('position_')

    @property
    def label(self) -> str:
        return self.__dict__.get('label_')

    def __setattr__(self, name, value):
        if name in self.__config__.alias_to_fields:
            name = self.__config__.alias_to_fields[name]
        else:
            warnings.warn(
                'Use of internal fields is discouraged. '
                'Please use external fields to set attributes.'
            )

        super().__setattr__(name, value)

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self is other

    def __repr__(self):
        return f'<{self.__class__.__name__}, id {id(self)}>'

    @validator('position_')
    def is_valid_position(cls, position):
        """Validator for attribute position"""
        if position is None:
            return None

        if not isinstance(position, u.unyt_array):
            warnings.warn('Positions are assumed to be in nm')
            position *= u.nm

        input_unit = position.units

        position = np.asarray(position, dtype=float, order='C')
        position = np.reshape(position, newshape=(3,), order='C')

        position *= input_unit
        position.convert_to_units(u.nm)
        return position

    @root_validator(pre=True)
    def inject_name(cls, values):
        if not values.get('name'):
            values['name'] = cls.__name__
        return values

    @classmethod
    def __new__(cls, *args: Any, **kwargs: Any) -> SiteT:
        if cls is Site:
            raise TypeError('Cannot instantiate abstract class of type Site')
        else:
            return object.__new__(cls)

    @classmethod
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__()
        apply_docs(cls, map_names=True)

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
