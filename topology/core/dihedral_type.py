import unyt as u

from topology.core.potential import Potential


class DihedralType(Potential):
    """A Potential between 4-bonded partners.

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
                 name='DihedralType',
                 expression='k * (1 + cos(n * phi - phi_eq))**2',
                 parameters={
                     'k': 1000 * u.Unit('kJ / (deg**2)'),
                     'theta_eq': 180 * u.deg,
                     'n': 1*u.dimensionless
                 },
                 independent_variables={'phi'}):

        super(DihedralType, self).__init__(name=name, expression=expression,
                parameters=parameters, independent_variables=independent_variables)

    def __repr__(self):
        return "<DihedralType {}, id {}>".format(self.name, id(self))

