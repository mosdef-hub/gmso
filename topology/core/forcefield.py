from topology.core.atom_type import AtomType
from topology.core.bond_type import BondType
from topology.core.angle_type import AngleType

from topology.exceptions import TopologyError

class Forcefield(object):
    """ A force field

    Parameters
    ---------
    name : str, optional, default 'Forcefield'
    """
    def __init__(self, name='Forcefield'):
        self.name = name
        self.atom_types = {}
        self.bond_types = {} 
        self.angle_types = {} 
        self.dihedral_types = {} 

    def __repr__(self):
        descr = list('<Forcefield ')
        descr.append(self.name + ' ')
        descr.append('{:d} AtomTypes, '.format(len(self.atom_types)))
        descr.append('{:d} BondTypes, '.format(len(self.bond_types)))
        descr.append('{:d} AngleTypes, '.format(len(self.angle_types)))
        descr.append('id: {}>'.format(id(self)))

        return ''.join(descr)

    def __setitem__(self, key, val):
        n_types = len(key.split('-'))
        if n_types == 1:
            if isinstance(val, AtomType):
                self.atom_types[key] = val
            else:
                raise TopologyError("Building Forcefield " + 
                        "non-AtomType parameters {}".format(val) +
                        " for AtomType {}".format(key))
        elif n_types == 2:
            if isinstance(val, BondType):
                self.bond_types[key] = val
            else:
                raise TopologyError("Building Forcefield " + 
                        "non-BondType parameters {}".format(val) +
                        " for BondType {}".format(key))
        elif n_types == 3:
            if isinstance(val, AngleType):
                self.angle_types[key] = val
            else:
                raise TopologyError("Building Forcefield " + 
                        "non-AngleType parameters {}".format(val) +
                        " for AngleType {}".format(key))
        else:
            raise TopologyError("Forcefield does not understand " + 
                    "parameter key {}".format(key))

    def __getitem__(self, key):
        n_types = len(key.split('-'))
        if n_types == 1:
            return self.atom_types.get(key, None)
        elif n_types == 2:
            return self.bond_types.get(key, None)
        elif n_types == 3:
            return self.angle_types.get(key, None)
        else:
            raise TopologyError("Forcefield does not understand " + 
                    "parameter key {}".format(key))

