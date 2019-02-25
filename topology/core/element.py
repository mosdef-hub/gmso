class Element(object):
    """An element."""
    def __init__(self, name=None, symbol=None, mass=None):
        self.name = name
        self.symbol = symbol
        self.mass = mass

    def lookup_from_name(self, name):
        return [elem for elem in elements if elem.name == name][0]

Hydrogen = Element(name='hydrogen', symbol='H', mass=1.007947)
Carbon = Element(name='carbon', symbol='C', mass=12.011)
Oxygen = Element(name='oxygen', symbol='O', mass=15.999)

elements = [
    Hydrogen,
    Carbon,
    Oxygen,
]
