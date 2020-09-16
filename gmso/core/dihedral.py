from typing import Tuple, Optional

from pydantic import Field

from gmso.abc.abstract_connection import Connection
from gmso.core.atom import Atom
from gmso.core.dihedral_type import DihedralType


class Dihedral(Connection):
    __base_doc__ = """A 4-partner connection between sites.

    This is a subclass of the gmso.Connection superclass.
    This class has strictly 4 members in its connection_members.
    The connection_type in this class corresponds to gmso.DihedralType.
    The connectivity of a dihedral is:

           m1–m2–m3–m4

    where m1, m2, m3, and m4 are connection members 1-4, respectively.

    Notes
    -----
    Inherits some methods from Connection:
        __eq__, __repr__, _validate methods
    Additional _validate methods are presented
    """
    connection_members_: Tuple[Atom, Atom, Atom, Atom] = Field(
        ...,
        description='The 4 atoms involved in the dihedral.'
    )

    dihedral_type_: Optional[DihedralType] = Field(
        default=None,
        description='DihedralType of this dihedral.'
    )

    @property
    def dihedral_type(self):
        return self.__dict__.get('dihedral_type_')

    @property
    def connection_type(self):
        # ToDo: Deprecate this?
        return self.__dict__.get('dihedral_type_')

    def __setattr__(self, key, value):
        if key == 'connection_type':
            super(Dihedral, self).__setattr__('dihedral_type', value)
        else:
            super(Dihedral, self).__setattr__(key, value)

    class Config:
        fields = {
            'dihedral_type_': 'dihedral_type',
            'connection_members_': 'connection_members'
        }
        alias_to_fields = {
            'dihedral_type': 'dihedral_type_',
            'connection_members': 'connection_members_'
        }
