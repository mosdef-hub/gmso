import warnings
import unyt as u

from gmso.core.potential import Potential
from gmso.exceptions import GMSOError
from gmso.utils._constants import IMPROPER_TYPE_DICT

class ImproperType(Potential):
    """A descripton of the interaction between 4 bonded partners.

    This is a subclass of the gmso.core.Potential superclass.

    ImproperType represents a improper type and includes the functional form
    describing its interactions. The functional form of the potential is stored
    as a `sympy` expression and the parameters, with units, are stored
    explicitly.  The AtomTypes that are used to define the improper type are
    stored as `member_types`.

    Parameters
    ----------
    name : str
    expression : str or sympy.Expression
        See `Potential` documentation for more information
    parameters : dict {str, unyt.unyt_quantity}
        See `Potential` documentation for more information
    independent vars : set of str
        See `Potential` documentation for more information
    member_types : list-like of str
        List-like of of gmso.AtomType.name defining the members of this
        improper type
    topology: gmso.core.Topology, default=None
        The topology of which this improper type is a part of
    set_ref: str
        The string name of the bookkeeping set in topology class.

    Notes
    ----
    Inherits many functions from gmso.Potential:
        __eq__, _validate functions

    """

    def __init__(self,
                 name='ImproperType',
                 expression='0.5 * k * ((phi - phi_eq))**2',
                 parameters=None,
                 independent_variables=None,
                 member_types=None,
                 topology=None,
                 set_ref='improper_type_set'):
        if parameters is None:
            parameters = {
                'k': 1000 * u.Unit('kJ / (deg**2)'),
                'phi_eq': 0 * u.deg,
            }
        if independent_variables is None:
            independent_variables = {'phi'}

        if member_types is None:
            member_types = list()

        super(ImproperType, self).__init__(name=name, expression=expression,
                                           parameters=parameters, independent_variables=independent_variables,
                                           topology=topology)
        self._set_ref = IMPROPER_TYPE_DICT
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
            warnings.warn("Changing an ImproperType's constituent "
                          "member types: {} to {}".format(self.member_types, val))
        self._member_types = _validate_four_member_type_names(val)

    def __repr__(self):
        return "<ImproperType {}, id {}>".format(self.name, id(self))


def _validate_four_member_type_names(types):
    """Ensure exactly 4 partners are involved in ImproperType"""
    if len(types) != 4 and len(types) != 0:
        raise GMSOError("Trying to create an ImproperType "
                            "with {} constituent types".format(len(types)))
    if not all([isinstance(t, str) for t in types]):
        raise GMSOError("Types passed to ImproperType "
                            "need to be strings corresponding to AtomType names")

    return types
