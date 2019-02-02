class Site(object):
    """A general site."""
    def __init__(self, name, position=None, element=None, atom_type=None):
        self.name = name
        if not position:
            self.position = np.zeros(3)
        else:
            self.position = position
        if element:
            self.element = element
        if atom_type:
            self.atom_type = atom_type
