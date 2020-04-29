import warnings
import unyt as u

from gmso.core.potential import ParametricPotential
from gmso.exceptions import GMSOError
from gmso.utils.decorators import confirm_dict_existence
from gmso.utils._constants import ANGLE_TYPE_DICT


class AngleType(ParametricPotential):
    """A descripton of the interaction between 3 bonded partners.

    This is a subclass of the gmso.core.Potential superclass.

    AngleType represents an angle type and includes the functional form
    describing its interactions. The functional form of the potential is stored
    as a `sympy` expression and the parameters, with units, are stored
    explicitly.  The AtomTypes that are used to define the angle type are
    stored as `member_types`.

    Parameters
    ----------
    name : str, optional
        A name uniquely representing the angle type
    expression : str or sympy.Expression
        See `Potential` documentation for more information
    parameters : dict {str, unyt.unyt_quantity}
        See `Potential` documentation for more information
    independent vars : set of str
        See `Potential` documentation for more information
    member_types : list-like of str
        List-like of of gmso.AtomType.name defining the members of this
        angle type
    gmso: gmso.core.Topology, default=None
        The Topology of which this angle_type is a part of
    set_ref: str
        The string name of the bookkeeping set in gmso class.

    Notes
    ----
    Inherits many functions from gmso.Potential:
        __eq__, _validate functions

    """

    def __init__(self,
                 name='AngleType',
                 expression='0.5 * k * (theta-theta_eq)**2',
                 parameters=None,
                 independent_variables=None,
                 member_types=None,
                 topology=None):
        if parameters is None:
            parameters = {
                'k': 1000 * u.Unit('kJ / (deg**2)'),
                'theta_eq': 180 * u.deg
            }
        if independent_variables is None:
            independent_variables = {'theta'}

        if member_types is None:
            member_types = list()
        super(AngleType, self).__init__(name=name, expression=expression,
                                        parameters=parameters, independent_variables=independent_variables,
                                        topology=topology)
        self._member_types = _validate_three_member_type_names(member_types)
        self._set_ref = ANGLE_TYPE_DICT

    @property
    def set_ref(self):
        return self._set_ref

    @property
    def member_types(self):
        return self._member_types

    @member_types.setter
    @confirm_dict_existence
    def member_types(self, val):
        if self.member_types != val:
            warnings.warn("Changing an AngleType's constituent "
                          "member types: {} to {}".format(self.member_types, val))
        self._member_types = _validate_three_member_type_names(val)

    def __repr__(self):
        return "<AngleType {}, id {}>".format(self.name, id(self))


def _validate_three_member_type_names(types):
    """Ensure exactly 3 partners are involved in AngleType"""
    if len(types) != 3 and len(types) != 0:
        raise GMSOError("Trying to create an AngleType "
                            "with {} constituent types".format(len(types)))
    if not all([isinstance(t, str) for t in types]):
        raise GMSOError("Types passed to AngleType "
                            "need to be strings corresponding to AtomType names")

    return types
