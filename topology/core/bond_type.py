import unyt as u
import sympy

from topology.core.potential import Potential
from topology.exceptions import TopologyError

class BondType(Potential):
    """A Potential between 2-bonded partners.
    
    Parameters
    ----------
    name : str
    expression : str or sympy.Expression
    parameters : dict
        {str, u.Unit}
    independent vars : set of str
        For canonical purposes, r needs to be an indpeendnet variable

    Notes
    ----
    Inherits many functions from topology.Potential:
        __eq__, _validate functions
    But we have specified an addiitonal validate method for canonicalizing
        the independent varaible
    """

    def __init__(self,
                 name='BondType',
                 expression='0.5 * k * (r-r_eq)**2',
                 parameters={
                     'k': 1000 * u.Unit('kJ / (nm**2)'),
                     'r_eq': 0.14 * u.nm
                 },
                 independent_variables={'r'}):

        super(BondType, self).__init__(name=name, expression=expression,
                parameters=parameters, independent_variables=independent_variables)

    def __repr__(self):
        return "<BondType {}, id {}>".format(self.name, id(self))

