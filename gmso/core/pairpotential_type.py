import unyt as u
import warnings

from gmso.core.potential import Potential
from gmso.utils.decorators import confirm_dict_existence
from gmso.exceptions import GMSOError
from gmso.utils._constants import BOND_TYPE_DICT


class PairPotentialType(Potential):
    """A descripton of the interaction between 2 bonded partners.

    This is a subclass of the gmso.core.Potential superclass.

    BondType represents a bond type and includes the functional form describing
    its interactions. The functional form of the potential is stored as a
    `sympy` expression and the parameters, with units, are stored explicitly.
    The AtomTypes that are used to define the bond type are stored as
    `member_types`.

    Parameters
    ----------
    name : str
        The name of the potential.
    expression : str or sympy.Expression
        See `Potential` documentation for more information
    parameters : dict {str, unyt.unyt_quantity}
        See `Potential` documentation for more information
    independent vars : set of str
        see `Potential` documentation for more information
    member_types : list-like of str
        List-like of of gmso.AtomType.name defining the members of this
        bond type

    Notes
    ----
    Inherits many functions from gmso.Potential:
        __eq__, _validate functions

    """

    def __init__(self,
                 name='PairPotentialType',
                 expression='4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
                 parameters=None,
                 independent_variables=None,
                 member_types=None,
                 topology=None,
                 set_ref='pairpotential_type_set'):
        if parameters is None:
            parameters = {'sigma': 0.3 * u.nm,
                          'epsilon': 0.3 * u.Unit('kJ')}
        if independent_variables is None:
            independent_variables = {'r'}

        if member_types is None:
            member_types = list()

        super(PairPotentialType, self).__init__(name=name, expression=expression,
                                       parameters=parameters, independent_variables=independent_variables,
                                       topology=topology)
        self._set_ref = PAIRPOTENTIAL_TYPE_DICT
        self._member_types = _validate_two_member_type_names(member_types)

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
            warnings.warn("Changing a PairPotentialType's constituent "
                          "member types: {} to {}".format(self.member_types, val))
        self._member_types = _validate_two_member_type_names(val)

    def __repr__(self):
        return "<PairPotentialType {}, id {}>".format(self.name, id(self))


def _validate_two_member_type_names(types):
    """Ensure exactly 2 partners are involved in PairPotentialType"""
    if len(types) != 2 and len(types) != 0:
        raise GMSOError("Trying to create a PairPotentialType"
                            "with {} constituent types".format(len(types)))
    if not all([isinstance(t, str) for t in types]):
        raise GMSOError("Types passed to PairPotentialType "
                            "need to be strings corresponding to AtomType names")

    return types
