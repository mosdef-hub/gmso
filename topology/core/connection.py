from topology.core.connection_type import ConnectionType

class Connection(object):
    """A simple connection between sites."""
    def __init__(self, site1=None, site2=None, update=True, connection_type=None):
        if site1:
            self.site1 = site1
        if site2:
            self.site2 = site2
        if update:
            site1.add_connection(site2)
            site2.add_connection(site1)

        if isinstance(connection_type, ConnectionType):
            self._connection_type = connection_type
        elif connection_type is None:
            self._connection_type = None
        else:
            self._connection_type = ConnectionType(connection_type)

    @property
    def connection_type(self):
        return self._connection_type

    @connection_type.setter
    def connection_type(self, connection_type):
        if isinstance(connection_type, ConnectionType):
            self._connection_type = connection_type
        elif connection_type is None:
            self._connection_type = None
        else:
            self._connection_type = ConnectionType(connection_type)

    def __eq__(self, other):
        # No comparison of connection_type in case of non-parametrization
        return ((self.site1 == other.site1 and self.site2 == other.site2) or
                (self.site2 == other.site1 and self.site1 == other.site2))
