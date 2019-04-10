import unyt as u

from topology.core.potential import Potential
from topology.exceptions import TopologyError


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
    member_types : list of topology.AtomType.name (str)

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
                 member_types=[]):

        super(BondType, self).__init__(name=name, expression=expression,
                parameters=parameters, independent_variables=independent_variables)

        self._member_types = _validate_two_member_type_names(member_types)

    @property
    def member_types(self):
        return self._member_types

    @member_types.setter
    def member_types(self, val):
       self._member_types = _validate_two_member_type_names(val)

    def __repr__(self):
        return "<BondType {}, id {}>".format(self.name, id(self))

def _validate_two_member_type_names(types):
    """Ensure 2 partners are involved in BondType"""
    if len(types) != 2 and len(types) != 0:
        raise TopologyError("Trying to create a BondType "
                "with {} constituent types". format(len(types)))
    if not all([isinstance(t, str) for t in types]):
        raise TopologyError("Types passed to BondType "
                            "need to be strings corresponding to AtomType names")

    return types

