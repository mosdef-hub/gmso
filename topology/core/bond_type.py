import unyt as u
import warnings

from topology.core.potential import Potential
from topology.utils.decorators import confirm_dict_existence
from topology.exceptions import TopologyError
from topology.utils._constants import BOND_TYPE_DICT


class BondType(Potential):
    """A Potential between 2-bonded partners.

    Parameters
    ----------
    name : str
    expression : str or sympy.Expression
        See `Potential` documentation for more information
    parameters : dict {str, unyt.unyt_quantity}
        See `Potential` documentation for more information
    independent vars : set of str
        see `Potential` documentation for more information
    member_types : list of topology.AtomType.name (str)

    Notes
    ----
    Inherits many functions from topology.Potential:
        __eq__, _validate functions
        """

    def __init__(self,
                 name='BondType',
                 expression='0.5 * k * (r-r_eq)**2',
                 parameters=None,
                 independent_variables=None,
                 member_types=None,
                 topology=None,
                 set_ref='bond_type_set'):
        if parameters is None:
            parameters = {
                'k': 1000 * u.Unit('kJ / (nm**2)'),
                'r_eq': 0.14 * u.nm
            }
        if independent_variables is None:
            independent_variables = {'r'}

        if member_types is None:
            member_types = list()

        super(BondType, self).__init__(name=name, expression=expression,
                                       parameters=parameters, independent_variables=independent_variables,
                                       topology=topology)
        self._set_ref = BOND_TYPE_DICT
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
            warnings.warn("Changing a BondType's constituent "
                          "member types: {} to {}".format(self.member_types, val))
        self._member_types = _validate_two_member_type_names(val)

    def __repr__(self):
        return "<BondType {}, id {}>".format(self.name, id(self))


def _validate_two_member_type_names(types):
    """Ensure 2 partners are involved in BondType"""
    if len(types) != 2 and len(types) != 0:
        raise TopologyError("Trying to create a BondType "
                            "with {} constituent types".format(len(types)))
    if not all([isinstance(t, str) for t in types]):
        raise TopologyError("Types passed to BondType "
                            "need to be strings corresponding to AtomType names")

    return types
