import unyt as u
import sympy

from topology.core.potential import Potential
from topology.exceptions import TopologyError

class AngleType(Potential):
    """A Potential between 3-bonded partners.
    
    Parameters
    ----------
    name : str
    expression : str or sympy.Expression
    parameters : dict
        {str, u.Unit}
    independent vars : set of str

    Notes
    ----
    Inherits many functions from topology.Potential:
        __eq__, _validate functions
    But we have specified an addiitonal validate method for canonicalizing
        the independent varaible
    """

    def __init__(self,
                 name='AngleType',
                 expression='0.5 * k * (theta-theta_eq)**2',
                 parameters={
                     'k': 1000 * u.Unit('kJ / (deg**2)'),
                     'r_eq': 0.14 * u.deg
                 },
                 independent_variables={'theta'}):

        super(AngleType, self).__init__(name=name, expression=expression,
                parameters=parameters, independent_variables=independent_variables)

    def __repr__(self):
        return "<AngleType {}, id {}>".format(self.name, id(self))

