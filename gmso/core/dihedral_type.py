import warnings
import unyt as u

from gmso.core.potential import Potential
from gmso.exceptions import GMSOError
from gmso.utils._constants import DIHEDRAL_TYPE_DICT

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
    member_types : list of gmso.AtomType.name (str)
    topology: gmso.core.Topology, the topology of which this dihedral type is a part of, default=None
    set_ref: (str), the string name of the bookkeeping set in topology class.

    Notes
    ----
    Inherits many functions from gmso.Potential:
        __eq__, _validate functions
    """

    def __init__(self,
                 name='DihedralType',
                 expression='k * (1 + cos(n * phi - phi_eq))**2',
                 parameters=None,
                 independent_variables=None,
                 member_types=None,
                 topology=None,
                 set_ref='dihedral_type_set'):
        if parameters is None:
            parameters = {
                'k': 1000 * u.Unit('kJ / (deg**2)'),
                'phi_eq': 180 * u.deg,
                'n': 1 * u.dimensionless
            }
        if independent_variables is None:
            independent_variables = {'phi'}

        if member_types is None:
            member_types = list()

        super(DihedralType, self).__init__(name=name, expression=expression,
                                           parameters=parameters, independent_variables=independent_variables,
                                           topology=topology)
        self._set_ref = DIHEDRAL_TYPE_DICT
        self._member_types = _validate_four_member_type_names(member_types)

    @property
    def set_ref(self):
        return self._set_ref

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
        raise GMSOError("Trying to create an DihedralType "
                            "with {} constituent types".format(len(types)))
    if not all([isinstance(t, str) for t in types]):
        raise GMSOError("Types passed to DihedralType "
                            "need to be strings corresponding to AtomType names")

    return types
