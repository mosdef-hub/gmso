import warnings
from gmso.core.potential import Potential
from gmso.core.site import Site
from gmso.exceptions import GMSOError


class Connection(object):
    """An abstract class that stores data about connections between sites.

    This class functions as a super-class for any connected groups (bonds, angles, dihedrals, etc).
    Each instance will have a property for the conection_type (bond_type, angle_type, dihedral_type)

    Parameters
    ----------

    connection_members : list of gmso.Site
        A list of constituents in this connection, in order.
    connection_type : gmso.Potential
        An instance of gmso.Potential that describes the potential, function and parameters of this interaction
    name : str, optional, default="Connection"
        A unique name for the connection. Used for writing hoomdxml bonds/angles/dihedrals.
        """
    def __init__(self, connection_members=None, connection_type=None, name="Connection"):
        if connection_members is None:
            connection_members = tuple()

        self._connection_members = _validate_connection_members(connection_members)
        self._connection_type = _validate_connection_type(connection_type)
        self._name = _validate_name(name)
        self._update_members()

    @property
    def connection_members(self):
        return self._connection_members

    @connection_members.setter
    def connection_members(self, connection_members):
        self._connection_members = _validate_connection_members(connection_members)

    @property
    def connection_type(self):
        return self._connection_type

    @connection_type.setter
    def connection_type(self, contype):
        self._connection_type = _validate_connection_type(contype)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, conname):
        self._name = _validate_name(conname)

    def equivalent_members(self):
        """Get a set of the equivalent connection member tuples

        Returns
        _______
        frozenset
            A unique set of tuples of equivalent connection members

        Notes
        _____
        For a bond:
            i, j == j, i
        where i and j are the connection members.
        """
        return frozenset(
                tuple(self.connection_members),
                tuple(reversed(self.connection_members))
                )

    def _equivalent_member_hash(self):
        """Get a unique hash representing the connection

        Returns
        _______
        int
            A unique hash to represent the connection members

        Notes
        _____
        Generalized for all connections, this is just a hashed tuple
        of the members. For specific connections (i.e. Bonds, Angles,
        Dihedrals, and Impropers), this function is overridden.
        """

        return hash(tuple(self.connection_members))

    def _update_members(self):
        for partner in self.connection_members:
            if self not in partner.connections:
                partner.add_connection(self)

    def __repr__(self):
        descr = '<{}-partner Connection, id {}, '.format(
                len(self.connection_members), id(self))
        descr += 'type {}'.format(self.connection_type)
        if self.name:
            descr += ', name {}'.format(self.name)
        descr += '>'

        return descr


def _validate_connection_members(connection_members):
    """Ensure all elements of entered connection_members are gmso.Site"""
    for partner in connection_members:
        if not isinstance(partner, Site):
            raise GMSOError("Supplied non-Site {}".format(partner))

    if len(set(connection_members)) != len(connection_members):
        raise GMSOError("Error, cannot add connection between same sites.")
    return tuple(connection_members)


def _validate_connection_type(c_type):
    """Ensure given connection_type is the gmso.Potential"""
    if c_type is None:
        warnings.warn("Non-parametrized Connection detected")
    elif not isinstance(c_type, Potential):
        raise GMSOError("Supplied non-Potential {}".format(c_type))
    return c_type


def _validate_name(conname):
    """Ensure given name is a string"""
    if not isinstance(conname, str):
        raise GMSOError("Supplied name {} is not a string".format(conname))
    return conname
