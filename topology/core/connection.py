import warnings
from topology.core.potential import Potential
from topology.core.site import Site
from topology.exceptions import TopologyError

class Connection(object):
    """ An abstract object that lists bonded partners and their type
    This functions as a super-class for any bonded groups (bonds,
    angles, dihedrals, etc), with a property for the conection_typie
    
    Parameters
    ----------
    connected_members : list of topology.Site
        A list of constituents in this bond. Should be in order
    connection_type : topology.Potential
        An instance of topology.Potential that describes
        the potential function and parameters of this interaction
        """
    def __init__(self, connected_members=[], connection_type=None):
        self._connected_members = _validate_connected_members(connected_members)
        self._connection_type = _validate_connection_type(connection_type)
        self._update_partners()

    @property
    def connected_members(self):
        return self._connected_members

    @connected_members.setter
    def connected_members(self, connected_members):
        self._connected_members = _validate_connected_members(connected_members)

    @property
    def connection_type(self):
        return self._connection_type

    @connection_type.setter
    def connection_type(self, contype):
        self_connection_type = _validate_connection_type(contype)

    def _update_partners(self):
        for partner in self.connected_members:
            if self not in partner.connections:
                partner.add_connection(self)

    def __repr__(self):
        descr = '<{}-partner Connection, id {}, '.format(
                len(self.connected_members), id(self))
        descr += 'type {}>'.format(self.connection_type)
        return descr

    def __eq__(self, other):
        bond_partner_match = (self.connected_members == other.connected_members)
        ctype_match = (self.connection_type == other.connection_type)
        return all([bond_partner_match, ctype_match])


def _validate_connected_members(connected_members):
    for partner in connected_members:
        if not isinstance(partner, Site):
            raise TopologyError("Supplied non-Site {}".format(partner))
    return connected_members

def _validate_connection_type(c_type):
    if c_type is None:
        warnings.warn("Non-parametrized Connection detected")
    elif not isinstance(c_type, Potential):
        raise TopologyError("Supplied non-Potential {}".format(c_type))
    return c_type
