"""Support for 3-partner connections between gmso.core.Atoms."""
from typing import Callable, ClassVar, Optional, Tuple

from pydantic import ConfigDict, Field

from gmso.abc.abstract_connection import Connection
from gmso.core.angle_type import AngleType
from gmso.core.atom import Atom


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
        ..., description="The 3 atoms involved in the angle."
    )

    angle_type_: Optional[AngleType] = Field(
        default=None, description="AngleType of this angle."
    )

    restraint_: Optional[dict] = Field(
        default=None,
        description="""
        Restraint for this angle, must be a dict with the following keys:
        'k' (unit of energy/mol), 'theta_eq' (unit of angle), 'n' (multiplicity, unitless).
        Refer to https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html
        for more information.
        """,
    )

    @property
    def angle_type(self):
        """Return the angle type if the angle is parametrized."""
        return self.__dict__.get("angle_type_")

    @property
    def connection_type(self):
        """Return the angle type if the angle is parametrized."""
        return self.__dict__.get("angle_type_")

    @property
    def restraint(self):
        """Return the restraint of this angle."""
        return self.__dict__.get("restraint_")

    def equivalent_members(self):
        """Return a set of the equivalent connection member tuples.

        Returns
        -------
        frozenset
            A unique set of tuples of equivalent connection members

        Notes
        -----
        For an angle:

            i, j, k == k, j, i

        where i, j and k are the connection members.
        """
        return frozenset(
            [self.connection_members, tuple(reversed(self.connection_members))]
        )

    def __setattr__(self, key, value):
        """Set the attributes of the angle."""
        if key == "connection_type":
            super(Angle, self).__setattr__("angle_type", value)
        else:
            super(Angle, self).__setattr__(key, value)

    # TODO[pydantic]: The following keys were removed: `fields`.
    # Check https://docs.pydantic.dev/dev-v2/migration/#changes-to-config for more information.
    model_config = ConfigDict(
        fields={
            "connection_members_": "connection_members",
            "angle_type_": "angle_type",
            "restraint_": "restraint",
        },
        alias_to_fields={
            "connection_members": "connection_members_",
            "angle_type": "angle_type_",
            "restraint": "restraint_",
        },
    )
