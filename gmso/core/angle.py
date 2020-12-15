from typing import Tuple, Optional, ClassVar, Callable

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
    __members_creator__: ClassVar[Callable] = Atom.parse_obj

    connection_members_: Tuple[Atom, Atom, Atom] = Field(
        ...,
        description='The 3 atoms involved in the angle.'
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

    def equivalent_members(self):
        """Get a set of the equivalent connection member tuples
        Returns
        _______
        frozenset
            A unique set of tuples of equivalent connection members
        Notes
        _____
        For an angle:
            i, j, k == k, j, i
        where i, j and k are the connection members.
        """

        return frozenset([
                self.connection_members,
                tuple(reversed(self.connection_members))
                ])

    def _equivalent_members_hash(self):
        """Returns a unique hash representing the connection
        Returns
        _______
        int
            A unique hash to represent the connection members
        Notes
        _____
        For an angle:
            i, j, k == k, j, i
        where i, j, and k are the connection members.
        Here, j is fixed and i and k are replaceable.
        """

        return hash(tuple([
            self.connection_members[1],
            frozenset([self.connection_members[0],
                       self.connection_members[2]])
            ]))

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
