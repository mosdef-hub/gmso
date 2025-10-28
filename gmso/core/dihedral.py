from typing import Callable, ClassVar, Optional, Tuple

from pydantic import ConfigDict, Field, model_validator

from gmso.abc.abstract_connection import Connection
from gmso.core.atom import Atom
from gmso.core.bond import Bond
from gmso.core.dihedral_type import DihedralType


class Dihedral(Connection):
    """A 4-partner connection between sites.

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

    Additional _validate methods are presented.
    """

    __members_creator__: ClassVar[Callable] = Atom.model_validate

    connectivity: ClassVar[Tuple[Tuple[int]]] = ((0, 1), (1, 2), (2, 3))

    connection_members_: Tuple[Atom, Atom, Atom, Atom] = Field(
        ...,
        description="The 4 atoms involved in the dihedral.",
        alias="connection_members",
    )

    bonds_: Tuple[Bond, Bond, Bond] = Field(
        default=None,
        description="""
        List of connection bonds.
        Ordered to align with connection_members, such that self.bonds_[0] is
        the bond between (self.connection_members[0], self.connection_members[1]).
        """,
        alias="bonds",
    )

    dihedral_type_: Optional[DihedralType] = Field(
        default=None,
        description="DihedralType of this dihedral.",
        alias="dihedral_type",
    )

    restraint_: Optional[dict] = Field(
        default=None,
        description="""
        Restraint for this dihedral, must be a dict with the following keys:
        'k' (unit of energy/(mol * angle**2)), 'phi_eq' (unit of angle), 'delta_phi' (unit of angle).
        Refer to https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html
        for more information.
        """,
        alias="restraint",
    )

    model_config = ConfigDict(
        alias_to_fields=dict(
            **Connection.model_config["alias_to_fields"],
            **{
                "dihedral_type": "dihedral_type_",
                "restraint": "restraint_",
            },
        )
    )

    @property
    def dihedral_type(self):
        return self.__dict__.get("dihedral_type_")

    @property
    def connection_type(self):
        # ToDo: Deprecate this?
        return self.__dict__.get("dihedral_type_")

    @property
    def restraint(self):
        """Return the restraint of this dihedral."""
        return self.__dict__.get("restraint_")

    @property
    def bonds(self):
        """Return the bonds that makeup this dihedral.

        Connectivity is ((0,1), (1,2), (2,3))
        """
        return self.__dict__.get("bonds_")

    @property
    def bonds_orders(self):
        """Return the bond_order strings of this dihedral."""
        return "".join([b.bond_order for b in self.bonds])

    def equivalent_members(self):
        """Get a set of the equivalent connection member tuples

        Returns
        _______
        frozenset
            A unique set of tuples of equivalent connection members

        Notes
        _____
        For a dihedral:

            i, j, k, l == l, k, j, i

        where i, j, k, and l are the connection members.
        """
        return frozenset(
            [self.connection_members, tuple(reversed(self.connection_members))]
        )

    def __setattr__(self, key, value):
        if key == "connection_type":
            super(Dihedral, self).__setattr__("dihedral_type_", value)
        else:
            super(Dihedral, self).__setattr__(key, value)

    @model_validator(mode="before")
    @classmethod
    def set_dependent_value_default(cls, data):
        """Automatically set bonds for this dihedral if connection_members is defined."""
        if "bonds" not in data and "connection_members" in data:
            atoms = data["connection_members"]
            data["bonds"] = (
                Bond(connection_members=(atoms[0], atoms[1])),
                Bond(connection_members=(atoms[1], atoms[2])),
                Bond(connection_members=(atoms[2], atoms[3])),
            )
        return data
