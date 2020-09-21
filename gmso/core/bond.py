from typing import Tuple, Optional

from pydantic import Field

from gmso.abc.abstract_connection import Connection
from gmso.core.atom import Atom
from gmso.core.bond_type import BondType


class Bond(Connection):
    __base_doc__ = """A 2-partner connection between sites.

    This is a subclass of the gmso.abc.Connection superclass.
    This class has strictly 2 members in its connection_members.
    The connection_type in this class corresponds to gmso.BondType.

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods.
    Additional _validate methods are presented.
    """
    connection_members_: Tuple[Atom, Atom] = Field(
        ...,
        description='The 2 atoms involved in the bond.'
    )

    bond_type_: Optional[BondType] = Field(
        default=None,
        description='BondType of this bond.'
    )

    @property
    def bond_type(self):
        return self.__dict__.get('bond_type_')

    @property
    def connection_type(self):
        # ToDo: Deprecate this?
        return self.__dict__.get('bond_type_')

    def __setattr__(self, key, value):
        if key == 'connection_type':
            super(Bond, self).__setattr__('bond_type', value)
        else:
            super(Bond, self).__setattr__(key, value)

    class Config:
        fields = {
            'bond_type_': 'bond_type',
            'connection_members_': 'connection_members'
        }
        alias_to_fields = {
            'bond_type': 'bond_type_',
            'connection_members': 'connection_members_'
        }
