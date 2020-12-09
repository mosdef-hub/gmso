from typing import Tuple, Optional
import unyt as u

from pydantic import Field

from gmso.core.parametric_potential import ParametricPotential
from gmso.utils._constants import BOND_TYPE_DICT


class BondType(ParametricPotential):
    __base_doc__ = """A descripton of the interaction between 2 bonded partners.

    This is a subclass of the gmso.core.Potential superclass.

    BondType represents a bond type and includes the functional form describing
    its interactions. The functional form of the potential is stored as a
    `sympy` expression and the parameters, with units, are stored explicitly.
    The AtomTypes that are used to define the bond type are stored as
    `member_types`.

    Notes
    ----
    Inherits many functions from gmso.ParametricPotential:
        __eq__, _validate functions
    """

    member_types_: Optional[Tuple[str, str]] = Field(
        None,
        description='List-like of of gmso.AtomType.name or gmso.AtomType.atomclass '
                    'defining the members of this bond type'
    )

    def __init__(self,
                 name='BondType',
                 expression=None,
                 parameters=None,
                 independent_variables=None,
                 potential_expression=None,
                 member_types=None,
                 topology=None):
        if potential_expression is None:
            if expression is None:
                expression = '0.5 * k * (r-r_eq)**2'

            if parameters is None:
                parameters = {
                    'k': 1000 * u.Unit('kJ / (nm**2)'),
                    'r_eq': 0.14 * u.nm
                }
            if independent_variables is None:
                independent_variables = {'r'}

        super(BondType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables,
            potential_expression=potential_expression,
            topology=topology,
            member_types=member_types,
            set_ref=BOND_TYPE_DICT
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
