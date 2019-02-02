class Site(object):
    """A general site."""
    def __init__(self, name, element=None, atom_type=None):
        self.name = name
        if element:
            self.element = element
        if atom_type:
            self.atom_type = atom_type
