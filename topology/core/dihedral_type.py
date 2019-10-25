import warnings
import unyt as u

from topology.core.potential import Potential
from topology.exceptions import TopologyError


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
    member_types : list of topology.AtomType.name (str)

    Notes
    ----
    Inherits many functions from topology.Potential:
        __eq__, _validate functions
    """

    def __init__(self,
                 name='DihedralType',
                 expression='k * (1 + cos(n * phi - phi_eq))**2',
                 parameters=None,
                 independent_variables=None,
                 member_types=None):
        if parameters is None:
            parameters = {
                     'k': 1000 * u.Unit('kJ / (deg**2)'),
                     'phi_eq': 180 * u.deg,
                     'n': 1*u.dimensionless
                 }
        if independent_variables is None:
            independent_variables = {'phi'}

        if member_types is None:
            member_types = list()

        super(DihedralType, self).__init__(name=name, expression=expression,
                parameters=parameters, independent_variables=independent_variables)

        self._member_types = _validate_four_member_type_names(member_types)

    @property
    def member_types(self):
        return self._member_types

    @member_types.setter
    def member_types(self, val):
        if self.member_types != val:
            warnings.warn("Changing an DihedralType's constituent "
                    "member types: {} to {}".format(self.member_types, val))
        self._member_types = _validate_four_member_type_names(val)

    def __repr__(self):
        return "<DihedralType {}, id {}>".format(self.name, id(self))

def _validate_four_member_type_names(types):
    """Ensure 4 partners are involved in DihedralType"""
    if len(types) != 4 and len(types) != 0:
        raise TopologyError("Trying to create an DihedralType "
                "with {} constituent types". format(len(types)))
    if not all([isinstance(t, str) for t in types]):
        raise TopologyError("Types passed to DihedralType "
                            "need to be strings corresponding to AtomType names")

    return types
