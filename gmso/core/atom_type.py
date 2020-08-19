import warnings
from typing import Optional, Set

import unyt as u
from pydantic import Field, validator

from gmso.core.parametric_potential import ParametricPotential
from gmso.utils.misc import unyt_to_hashable
from gmso.utils.decorators import confirm_dict_existence
from gmso.utils._constants import ATOM_TYPE_DICT


class AtomType(ParametricPotential):
    __base_doc__ = """A description of non-bonded interacitons between sites.

    This is a subclass of the gmso.core.Potential superclass.

    AtomType represents an atom type and includes the functional form
    describing its interactions and, optionally, other properties such as mass
    and charge.  This class inhereits from Potential, which stores the
    non-bonded interaction between atoms or sites. The functional form of the
    potential is stored as a `sympy` expression and the parameters, with units,
    are stored explicitly.
    """

    mass_: Optional[u.unyt_array] = Field(
        0.0 * u.gram / u.mol,
        description='The mass of the atom type'
    )

    charge_: Optional[u.unyt_array] = Field(
        0.0 * u.elementary_charge,
        description='The charge of the atom type'
    )

    atomclass_: Optional[str] = Field(
        '',
        description='The class of the atomtype'
    )

    doi_: Optional[str] = Field(
        '',
        description='Digital Object Identifier of publication where this atom type was'
    )

    overrides_: Optional[Set[str]] = Field(
        set(),
        description=' Set of other atom types that this atom type overrides'
    )

    definition_: Optional[str] = Field(
        '',
        description='SMARTS string defining this atom type'
    )

    description_: Optional[str] = Field(
        '',
        description='Description for the AtomType'
    )

    def __init__(self,
                 name='AtomType',
                 mass=0.0 * u.gram / u.mol,
                 charge=0.0 * u.elementary_charge,
                 expression='4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
                 parameters=None,
                 independent_variables=None,
                 atomclass='', doi='', overrides=None, definition='',
                 description='',
                 topology=None):
        if parameters is None:
            parameters = {'sigma': 0.3 * u.nm,
                          'epsilon': 0.3 * u.Unit('kJ')}
        if independent_variables is None:
            independent_variables = {'r'}

        if overrides is None:
            overrides = set()

        kwargs = {
            'mass': mass,
            'charge': charge,
            'atomclass': atomclass,
            'doi': doi,
            'overrides': overrides,
            'description': description,
            'definition': definition,
            'set_ref': ATOM_TYPE_DICT
        }

        super(AtomType, self).__init__(
            name=name,
            expression=expression,
            parameters=parameters,
            independent_variables=independent_variables,
            topology=topology,
            **kwargs
        )

    @property
    def charge(self):
        return self.__dict__.get('charge_')

    @property
    def mass(self):
        return self.__dict__.get('mass_')

    @property
    def atomclass(self):
        return self.__dict__.get('atomclass_')

    @property
    def doi(self):
        return self.__dict__.get('doi_')

    @property
    def overrides(self):
        return self.__dict__.get('overrides_')

    @property
    def description(self):
        return self.__dict__.get('description_')

    @property
    def definition(self):
        return self.__dict__.get('definition_')

    def __hash__(self):
        return hash(
            tuple(
                (
                    self.name,
                    unyt_to_hashable(self.mass),
                    unyt_to_hashable(self.charge),
                    self.expression,
                    tuple(self.independent_variables),
                    tuple(self.parameters.keys()),
                    tuple(unyt_to_hashable(val) for val in self.parameters.values())
                )
            )
        )

    @validator('mass_', pre=True)
    def validate_mass(cls, mass):
        """Check to see that a mass is a unyt array of the right dimension"""
        if not isinstance(mass, u.unyt_array):
            warnings.warn("Masses are assumed to be g/mol")
            mass *= u.gram / u.mol
        elif mass.units.dimensions != (u.gram / u.mol).units.dimensions:
            warnings.warn("Masses are assumed to be g/mol")
            mass = mass.value * u.gram / u.mol
        else:
            pass

        return mass

    @validator('charge_', pre=True)
    def validate_charge(cls, charge):
        """Check to see that a charge is a unyt array of the right dimension"""
        if not isinstance(charge, u.unyt_array):
            warnings.warn("Charges are assumed to be elementary charge")
            charge *= u.elementary_charge
        elif charge.units.dimensions != u.elementary_charge.units.dimensions:
            warnings.warn("Charges are assumed to be elementary charge")
            charge = charge.value * u.elementary_charge
        else:
            pass

        return charge

    class Config:
        fields = {
            'mass_': 'mass',
            'charge_': 'charge',
            'atomclass_': 'atomclass',
            'overrides_': 'overrides',
            'doi_': 'doi',
            'description_': 'description',
            'definition_': 'definition',
        }

        alias_to_fields = {
            'mass': 'mass_',
            'charge': 'charge_',
            'atomclass': 'atomclass_',
            'overrides': 'overrides_',
            'doi': 'doi_',
            'description': 'description_',
            'definition': 'definition_',
        }
