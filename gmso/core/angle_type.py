from typing import Tuple, Optional
import unyt as u

from pydantic import Field

from gmso.core.parametric_potential import ParametricPotential
from gmso.utils._constants import ANGLE_TYPE_DICT


class AngleType(ParametricPotential):
    __base_doc__ = """A descripton of the interaction between 3 bonded partners.

    This is a subclass of the gmso.core.Potential superclass.

    AngleType represents an angle type and includes the functional form
    describing its interactions. The functional form of the potential is stored
    as a `sympy` expression and the parameters, with units, are stored
    explicitly.  The AtomTypes that are used to define the angle type are
    stored as `member_types`.

    Notes
    ----
    Inherits many functions from gmso.Potential:
        __eq__, _validate functions
    """

    member_types_: Optional[Tuple[str, str, str]] = Field(
        None,
        description='List-like of gmso.AtomType.name defining the members of this angle type'
    )

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

        super(AngleType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables,
            topology=topology,
            member_types=member_types,
            set_ref=ANGLE_TYPE_DICT
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
