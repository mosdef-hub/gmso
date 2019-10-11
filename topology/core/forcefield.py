from topology.core.atom_type import AtomType
from topology.core.bond_type import BondType
from topology.core.angle_type import AngleType

from topology.exceptions import TopologyError

class Forcefield(object):
    """ A force field

    Parameters
    ---------
    name : str, optional, default 'Forcefield'

    Notes
    -----
    While each set of parameters can be accessed via ff.atom_types,
    ff.bond_types, etc., they can also be accessed via
    ff['atomtype'], ff['bondtype'] for convenience and safer lookups

    Attributes
    ---------
    atom_type_definitions : dict
        'string':'AtomType'
    atom_types : dict
        'A':topology.AtomType
    bond_types : dict
        'A-B':topology.BondType
    angle_types : dict
        'A-B-C':topology.AngleType
    """
    def __init__(self, name='Forcefield'):
        self.name = name
        self.atom_type_definitions = {}
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
        types = key.split('-')
        if len(types) == 1:
            to_return = None
            return self.atom_types.get(key, to_return)
        elif len(types) == 2:
            to_return = None
            permutations = ['{first}-{second}'.format(first=types[0],
                                                        second=types[1]),
                            '{second}-{first}'.format(first=types[0],
                                                        second=types[1])]
            for permutation in permutations:
                to_return = self.bond_types.get(permutation, to_return)

            return to_return
        elif len(types)  == 3:
            to_return = None
            slots = {'first':types[0], 'second':types[1],
                    'third':types[2]}
            permutations  = ['{first}-{second}-{third}'.format(**slots),
                            '{third}-{second}-{first}'.format(**slots)]
            for permutation in permutations:
                to_return = self.angle_types.get(permutation, to_return)

            return to_return
        else:
            raise TopologyError("Forcefield does not understand " + 
                    "parameter key {}".format(key))

