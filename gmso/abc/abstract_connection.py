from typing import Sequence, Optional

from pydantic import Field, root_validator

from gmso.abc.abstract_site import Site
from gmso.abc.gmso_base import GMSOBase
from gmso.exceptions import GMSOError


class Connection(GMSOBase):
    __base_doc__ = """An abstract class that stores data about connections between sites.

    This class functions as a super-class for any connected groups (bonds, angles, dihedrals, etc).
    Each instance will have a property for the conection_type (bond_type, angle_type, dihedral_type)
    """

    name_: str = Field(
        default='',
        description='Name of the connection. Defaults to class name'
    )

    connection_members_: Optional[Sequence[Site]] = Field(
        default=None,
        description='A list of constituents in this connection, in order.'
    )

    @property
    def connection_members(self):
        return self.__dict__.get('connection_members_')

    @property
    def name(self):
        return self.__dict__.get('name_')

    @root_validator(pre=True)
    def validate_fields(cls, values):
        connection_members = values.get('connection_members')
        if not all(isinstance(x, Site) for x in connection_members):
            try:
                connection_members = [cls.__members_creator__(x) for x in connection_members]
            except:
                raise TypeError(
                    f'A non-site object provided to be a connection member'
                )
        if len(set(connection_members)) != len(connection_members):
            raise GMSOError(
                f'Trying to create a {cls.__name__} between '
                f'same sites. A {cls.__name__} between same '
                f'{type(connection_members[0]).__name__}s is not allowed'
            )
        if not values.get('name'):
            values['name'] = cls.__name__
        return values

    def __repr__(self):
        descr = '<{}-partner Connection, id {}, '.format(
                len(self.connection_members), id(self))
        descr += 'type {}'.format(self.connection_type)
        if self.name:
            descr += ', name {}'.format(self.name)
        descr += '>'

        return descr

    class Config:
        fields = {
            'name_': 'name',
            'connection_members_': 'connection_members'
        }

        alias_to_fields = {
            'name': 'name_',
            'connection_members': 'connection_members_'
        }
