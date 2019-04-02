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
    types : list of topology.AtomType.name (str)

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
                 independent_variables={'theta'},
                 types=[]):

        super(AngleType, self).__init__(name=name, expression=expression,
                parameters=parameters, independent_variables=independent_variables)

        self._types = _validate_three_atomtypes(types)

    @property
    def types(self):
        return self._types

    @types.setter
    def types(self, val):
       self._types = _validate_three_atomtypes(val)

    def __repr__(self):
        return "<AngleType {}, id {}>".format(self.name, id(self))

def _validate_three_atomtypes(types):
    """Ensure 3 partners are involved in BondType"""
    if len(types) != 3:
        raise TopologyError("Trying to create an AngleType "
                "with {} constituent types". format(len(types)))
    if not all([isinstance(t, str) for t in types]):
        raise TopologyError("Types passed to AngleType "
                            "need to be strings corresponding to AtomType names")

    return types

