from typing import Tuple, Optional

from pydantic import Field

from gmso.core.atom import Atom
from gmso.abc.abstract_connection import Connection
from gmso.core.angle_type import AngleType


class Angle(Connection):
    __base_doc__ = """A 3-partner connection between Atoms.

    This is a subclass of the gmso.Connection superclass.
    This class has strictly 3 members in its connection members.
    The connection_type in this class corresponds to gmso.AngleType.
    
    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods
    Additional _validate methods are presented
    """

    connection_members_: Tuple[Atom, Atom, Atom] = Field(
        ...,
        description='3 Atoms of an angle.'
    )

    angle_type_: Optional[AngleType] = Field(
        default=None,
        description='AngleType of this angle.'
    )

    @property
    def angle_type(self):
        return self.__dict__.get('angle_type_')

    @property
    def connection_type(self):
        return self.__dict__.get('angle_type_')

    def __setattr__(self, key, value):
        if key == 'connection_type':
            super(Angle, self).__setattr__('angle_type', value)
        else:
            super(Angle, self).__setattr__(key, value)

    class Config:
        fields = {
            'connection_members_': 'connection_members',
            'angle_type_': 'angle_type'
        }
        alias_to_fields = {
            'connection_members': 'connection_members_',
            'angle_type': 'angle_type_'
        }
