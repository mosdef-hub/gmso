from typing import Tuple, Optional

from pydantic import Field

from gmso.abc.abstract_connection import Connection
from gmso.core.atom import Atom
from gmso.core.improper_type import ImproperType


class Improper(Connection):
    __base_doc__ = """sA 4-partner connection between sites.

    This is a subclass of the gmso.Connection superclass.
    This class has strictly 4 members in its connection_members.
    The connection_type in this class corresponds to gmso.ImproperType
    The connectivity of an improper is:

         m2
         |
         m1
        /  \
        m3  m4

    where m1, m2, m3, and m4 are connection members 1-4, respectively.

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods
            
    Additional _validate methods are presented
    """
    connection_members_: Tuple[Atom, Atom, Atom, Atom] = Field(
        ...,
        description='The 4 atoms of this improper. Central site first, '
                    'then the three atoms connected to the central site.'
    )

    improper_type_: Optional[ImproperType] = Field(
        default=None,
        description='ImproperType of this improper.'
    )

    @property
    def improper_type(self):
        return self.__dict__.get('improper_type_')

    @property
    def connection_type(self):
        # ToDo: Deprecate this?
        return self.__dict__.get('improper_type_')

    def equivalent_members(self):
        """Get a set of the equivalent connection member tuples

        Returns
        -------
        frozenset
            A unique set of tuples of equivalent connection members

        Notes
        -----
        For an improper:

            i, j, k, l == i, k, j, l

        where i, j, k, and l are the connection members.
        """
        equiv_members = [self.connection_members[0],
                         self.connection_members[2],
                         self.connection_members[1],
                         self.connection_members[3]]

        return frozenset([
                self.connection_members,
                tuple(equiv_members)
                ])

    def _equivalent_members_hash(self):
        """Returns a unique hash representing the connection

        Returns
        -------
        int
            A unique hash to represent the connection members

        Notes
        -----
        For an improper:
            i, j, k, l == i, k, j, l
        where i, j, k, and l are the connection members.
        Here j and k are interchangeable and i and l are fixed.
        """

        return hash(tuple([
            self.connection_members[0],
            self.connection_members[3],
            frozenset([self.connection_members[1],
                       self.connection_members[2]])
            ]))


    def __setattr__(self, key, value):
        if key == 'connection_type':
            super(Improper, self).__setattr__('improper_type', value)
        else:
            super(Improper, self).__setattr__(key, value)

    class Config:
        fields = {
            'improper_type_': 'improper_type',
            'connection_members_': 'connection_members'
        }
        alias_to_fields = {
            'improper_type': 'improper_type_',
            'connection_members': 'connection_members_'
        }
