import numpy as np
import sympy
import unyt as u

class Connection(object):
    """ An abstract object that lists bonded partners and their type """
    def __init__(self, bond_partners=[], connection_type=None):
        self._bond_partners = _validate_bond_partners(bond_partners)
        self._connection_type = _validate_connection_type(connection_type)

    def __repr__(self):
        descr = '<{}-partner connection between'.format(len(self.bond_partners))
        descr += ' '.join([site for site in self.bond_partners])
        descr += ', type {}'.format(self.connection_type)
        return(descr)
