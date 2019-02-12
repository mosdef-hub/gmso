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

        self._connection_type = connection_type

    @property
    def connection_type(self):
        return self._connection_type

    @connection_type.setter
    def connection_type(self, val):
        self._connection_type = val

    def __eq__(self, other):
        # No comparison of connection_type in case of non-parametrization
        return ((self.site1 == other.site1 and self.site2 == other.site2) or
                (self.site2 == other.site1 and self.site1 == other.site2))
