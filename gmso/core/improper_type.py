from typing import Tuple, Optional
import unyt as u

from pydantic import Field

from gmso.core.parametric_potential import ParametricPotential
from gmso.utils._constants import IMPROPER_TYPE_DICT


class ImproperType(ParametricPotential):
    """A descripton of the interaction between 4 bonded partners.

    This is a subclass of the gmso.core.Potential superclass.

    ImproperType represents a improper type and includes the functional form
    describing its interactions. The functional form of the potential is stored
    as a `sympy` expression and the parameters, with units, are stored
    explicitly.  The AtomTypes that are used to define the improper type are
    stored as `member_types`.

    The connectivity of an improper is:

                   m2
                   |
                   m1
                  / \
                 m3  m4

    where m1, m2, m3, and m4 are connection members 1-4, respectively.

    Notes
    ----
    Inherits many functions from gmso.Potential:
        __eq__, _validate functions
    """

    member_types_: Optional[Tuple[str, str, str, str]] = Field(
        None,
        description='List-like of of gmso.AtomType.name or gmso.AtomType.atomclass '
                    'defining the members of this improper type'
    )

    def __init__(self,
                 name='ImproperType',
                 expression='0.5 * k * ((phi - phi_eq))**2',
                 parameters=None,
                 independent_variables=None,
                 member_types=None,
                 topology=None):
        if parameters is None:
            parameters = {
                'k': 1000 * u.Unit('kJ / (deg**2)'),
                'phi_eq': 0 * u.deg,
            }
        if independent_variables is None:
            independent_variables = {'phi'}

        super(ImproperType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables,
            topology=topology,
            member_types=member_types,
            set_ref=IMPROPER_TYPE_DICT
        )

    @property
    def member_types(self):
        return self.__dict__.get('member_types_')

    class Config:
        fields = {
            'member_types_': 'member_types'
        }

        alias_to_fields = {
            'member_types': 'member_types_'
        }
