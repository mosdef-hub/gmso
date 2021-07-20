"""Support for improper style connections (4-member connection)."""
from typing import Callable, ClassVar, Iterable, List, Optional, Tuple

from boltons.setutils import IndexedSet
from pydantic import Field, ValidationError, validator

from gmso.abc.abstract_connection import Connection
from gmso.core.atom import Atom
from gmso.core.improper_type import ImproperType
from gmso.utils.misc import validate_type


class BaseImproper(Connection):
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
    __members_creator__: ClassVar[Callable] = Atom.parse_obj

    connection_members_: Tuple[Atom, Atom, Atom, Atom] = Field(
        ...,
        description="The 4 atoms of this improper. Central atom first, "
        "then the three atoms connected to the central site.",
    )

    def equivalent_members(self):
        """Get a set of the equivalent connection member tuples.

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
        equiv_members = [
            self.connection_members[0],
            self.connection_members[2],
            self.connection_members[1],
            self.connection_members[3],
        ]

        return frozenset([self.connection_members, tuple(equiv_members)])

    def _equivalent_members_hash(self):
        """Return a unique hash representing the connection.

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
        return hash(
            tuple(
                [
                    self.connection_members[0],
                    self.connection_members[3],
                    frozenset(
                        [self.connection_members[1], self.connection_members[2]]
                    ),
                ]
            )
        )

    class Config:
        """Pydantic configuration to link fields to their public attribute."""

        fields = {
            "connection_members_": "connection_members",
        }
        alias_to_fields = {
            "connection_members": "connection_members_",
        }


class Improper(BaseImproper):
    improper_type_: Optional[ImproperType] = Field(
        default=None, description="ImproperType of this improper."
    )

    @property
    def improper_type(self):
        """Return Potential object for this connection if it exists."""
        return self.__dict__.get("improper_type_")

    @property
    def connection_type(self):
        """Return Potential object for this connection if it exists."""
        # ToDo: Deprecate this?
        return self.__dict__.get("improper_type_")

    def __setattr__(self, key, value):
        """Set attribute override to support connection_type key."""
        if key == "connection_type":
            super().__setattr__("improper_type", value)
        else:
            super().__setattr__(key, value)

    class Config:
        """Pydantic configuration to link fields to their public attribute."""

        fields = {
            "improper_type_": "improper_type",
        }
        alias_to_fields = {
            "improper_type": "improper_type_",
        }


class LayeredImproper(BaseImproper):
    improper_types_: Optional[List[Improper]] = Field(
        default=None, description="ImproperTypes of this improper."
    )

    @property
    def improper_types(self):
        return self.__dict__.get("improper_types_")

    @property
    def connection_types(self):
        # ToDo: Deprecate this?
        return self.__dict__.get("improper_types_")

    def __setattr__(self, key, value):
        if key == "connection_types":
            super().__setattr__("improper_types_", value)
        else:
            super().__setattr__(key, value)

    @validator("improper_types_", pre=True, always=True)
    def validate_dihedral_types(cls, improper_types):
        if not isinstance(improper_types, Iterable):
            raise ValidationError("ImproperTypes should be iterable", cls)

        validate_type(improper_types, ImproperType)
        return IndexedSet(improper_types)

    class Config:
        fields = {
            "improper_types_": "improper_types",
        }
        alias_to_fields = {
            "improper_types": "improper_types_",
        }
