import unyt as u

from topology.core.potential import Potential

class AngleType(Potential):
    """A Potential between 3-bonded partners.

    Parameters
    ----------
    name : str
    expression : str or sympy.Expression
        See `Potential` documentation for more information
    parameters : dict {str, unyt.unyt_quantity}
        See `Potential` documentation for more information
    independent vars : set of str
        See `Potential` documentation for more information

    Notes
    ----
    Inherits many functions from topology.Potential:
        __eq__, _validate functions
    """

    def __init__(self,
                 name='AngleType',
                 expression='0.5 * k * (theta-theta_eq)**2',
                 parameters={
                     'k': 1000 * u.Unit('kJ / (deg**2)'),
                     'theta_eq': 180 * u.deg
                 },
                 independent_variables={'theta'}):

        super(AngleType, self).__init__(name=name, expression=expression,
                parameters=parameters, independent_variables=independent_variables)

    def __repr__(self):
        return "<AngleType {}, id {}>".format(self.name, id(self))

