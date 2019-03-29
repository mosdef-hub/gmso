import unyt as u

from topology.core.potential import Potential
from topology.core.atom_type import AtomType


class BondType(Potential):
    """A Potential between 2-bonded partners.

    Parameters
    ----------
    name : str
    expression : str or sympy.Expression
        See `Potential` documentation for more information
    parameters : dict {str, unyt.unyt_quantity}
        See `Potential` documentation for more information
    independent vars : set of str
        see `Potential` documentation for more information
    types : list of topology.AtomType

    Notes
    ----
    Inherits many functions from topology.Potential:
        __eq__, _validate functions
        """

    def __init__(self,
                 name='BondType',
                 expression='0.5 * k * (r-r_eq)**2',
                 parameters={
                     'k': 1000 * u.Unit('kJ / (nm**2)'),
                     'r_eq': 0.14 * u.nm
                 },
                 independent_variables={'r'},
                 types=[]):

        super(BondType, self).__init__(name=name, expression=expression,
                parameters=parameters, independent_variables=independent_variables)

        self._types = _validate_two_atomtypes(types)

    @property
    def types(self):
        return self._types

    @types.setter
    def types(self, val):
       self._types = _validate_two_atomtypes(val)

    def __repr__(self):
        return "<BondType {}, id {}>".format(self.name, id(self))

def _validate_two_atomtypes(types):
    """Ensure 2 partners are involved in BondType"""
    if len(types) != 2:
        raise TopologyError("Trying to create a BondType "
                "with {} constituent types". format(len(types)))
    if not all([isinstance(t, AtomType) for t in types]):
        raise TopologyError("Types passed to BondType "
                            "need to be topology.AtomTypes")

    return types

