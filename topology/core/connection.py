class Connection(object):
    """A simple connection between sites."""
    def __init__(self, site1=None, site2=None):
        if site1:
            self.site1 = site1
        if site2:
            self.site2 = site2
        site1.add_connection(site2)
        site2.add_connection(site1)

