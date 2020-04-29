import warnings
import unyt as u

from gmso.core.potential import ParametricPotential
from gmso.exceptions import GMSOError
from gmso.utils._constants import DIHEDRAL_TYPE_DICT

class DihedralType(ParametricPotential):
    """A descripton of the interaction between 4 bonded partners.

    This is a subclass of the gmso.core.Potential superclass.

    DihedralType represents a dihedral type and includes the functional form
    describing its interactions. The functional form of the potential is stored
    as a `sympy` expression and the parameters, with units, are stored
    explicitly.  The AtomTypes that are used to define the dihedral type are
    stored as `member_types`.

    The connectivity of a dihedral is:

           m1–m2–m3–m4

    where m1, m2, m3, and m4 are connection members 1-4, respectively.

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
        dihedral type
    topology: gmso.core.Topology, default=None
        The topology of which this dihedral type is a part of
    set_ref: str
        The string name of the bookkeeping set in topology class.

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
                 topology=None):
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
                                           parameters=parameters,
                                           independent_variables=independent_variables,
                                           topology=topology,
                                           dict_ref=DIHEDRAL_TYPE_DICT)
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
    """Ensure exactly 4 partners are involved in DihedralType"""
    if len(types) != 4 and len(types) != 0:
        raise GMSOError("Trying to create an DihedralType "
                            "with {} constituent types".format(len(types)))
    if not all([isinstance(t, str) for t in types]):
        raise GMSOError("Types passed to DihedralType "
                            "need to be strings corresponding to AtomType names")

    return types
