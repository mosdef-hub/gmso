import warnings
import unyt as u

from topology.core.potential import Potential
from topology.exceptions import TopologyError
from topology.utils._constants import DIHEDRAL_TYPE_DICT

class VirtualSiteType(Potential):
    """A virtual site construction.

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
    topology: topology.core.Topology, the topology of which this virtual site type is a part of, default=None
    set_ref: (str), the string name of the bookkeeping set in topology class.

    Notes
    ----
    Inherits many functions from topology.Potential:
        __eq__, _validate functions
    """

    def __init__(self,
                 name='VirtualSiteType',
                 expression='1 - a',
                 parameters=None,
                 independent_variables=None,
                 member_types=None,
                 topology=None,
                 set_ref='virtualsite_type_set'):
        if parameters is None:
            parameters = {
                'a': 0.35 * u.dimensionless
            }

        if member_types is None:
            member_types = list()

        super(VirtualSiteType, self).__init__(name=name, expression=expression,
                                           parameters=parameters, independent_variables=independent_variables,
                                           topology=topology)
        #self._set_ref = VIRTUALSITE_TYPE_DICT
        self._member_types = member_types

    @property
    def set_ref(self):
        return self._set_ref

    @property
    def member_types(self):
        return self._member_types

    @member_types.setter
    def member_types(self, val):
        if self.member_types != val:
            warnings.warn("Changing an VirtualSiteType's constituent "
                          "member types: {} to {}".format(self.member_types, val))
        self._member_types = _validate_four_member_type_names(val)

    def __repr__(self):
        return "<VirtualSiteType {}, id {}>".format(self.name, id(self))
