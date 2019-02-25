class Element(object):
    """An element."""
    def __init__(self, name=None, symbol=None, mass=None):
        self.name = name
        self.symbol = symbol
        self.mass = mass

    def lookup_from_name(self, name):
        for elem in elements:
            if elem.name == name:
                return elem

    def __repr__(self):
        descr = list('<Element ')
        descr.append(self.name + ', ')
        descr.append(self.symbol + ', ')
        descr.append(str(self.mass) + ' amu>')
        return ''.join(descr)

Hydrogen = Element(name='hydrogen', symbol='H', mass=1.007947)
Carbon = Element(name='carbon', symbol='C', mass=12.011)
Oxygen = Element(name='oxygen', symbol='O', mass=15.999)

elements = [
    Hydrogen,
    Carbon,
    Oxygen,
]
