import warnings
from typing import Optional, Union

import unyt as u

from pydantic import Field, validator

from gmso.core.element import Element
from gmso.abc.abstract_site import Site
from gmso.core.atom_type import AtomType
from gmso.utils.misc import ensure_valid_dimensions
from gmso.utils._constants import UNIT_WARNING_STRING


class Atom(Site):
    __base_doc__ = """An atom represents a single element association in a topology.
    
    Atoms are the representation of an element within `gmso` that describes any general 
    atom in a molecular simulation. Atoms also contain information that are unique to 
    elements vs other types of interaction sites in molecular simulations. 
    For example, charge, mass, and periodic table information.

    Notes
    -----
    Atoms have all the attributes inherited from the base Site class,
    The order of precedence when attaining properties `charge` and `mass` is:
        
        1. atom.charge > atom.atom_type.charge 
        2. atom.mass > atom.atom_type.mass
        
    Examples
    --------
        >>> from gmso.core.atom import Atom
        >>> atom1 = Atom(name='lithium')

    See Also
    --------
    gmso.abc.AbstractSite
        An Abstract Base class for implementing site objects in GMSO. The class Atom bases from
        the gmso.abc.abstract site class
    """
    charge_: Optional[Union[u.unyt_quantity, float]] = Field(
        None,
        description='Charge of the atom',
    )

    mass_: Optional[Union[u.unyt_quantity, float]] = Field(
        None,
        description='Mass of the atom'
    )

    element_: Optional[Element] = Field(
        None,
        description='Element associated with the atom'
    )

    atom_type_: Optional[AtomType] = Field(
        None,
        description='AtomType associated with the atom'
    )

    @property
    def charge(self) -> Union[u.unyt_quantity, None]:
        charge = self.__dict__.get('charge_', None)
        atom_type = self.__dict__.get('atom_type_', None)
        if charge is not None:
            return charge
        elif atom_type is not None:
            return atom_type.charge
        else:
            return None

    @property
    def mass(self) -> Union[u.unyt_quantity, None]:
        mass = self.__dict__.get('mass_', None)
        atom_type = self.__dict__.get('atom_type_', None)
        if mass is not None:
            return mass
        elif atom_type is not None:
            return atom_type.mass
        else:
            return None

    @property
    def element(self) -> Union[Element, None]:
        return self.__dict__.get('element_', None)

    @property
    def atom_type(self) -> Union[AtomType, None]:
        return self.__dict__.get('atom_type_', None)

    def __le__(self, other):
        if isinstance(other, Atom):
            return hash(self) <= hash(other)
        else:
            raise TypeError(
                f'Cannot compare equality between {type(self)} and {type(other)}'
            )

    def __lt__(self, other):
        if isinstance(other, Atom):
            return hash(self) < hash(other)
        else:
            raise TypeError(
                f'Cannot compare equality between {type(self)} and {type(other)}'
            )

    @validator('charge_')
    def is_valid_charge(cls, charge):
        if charge is None:
            return None
        if not isinstance(charge, u.unyt_array):
            warnings.warn(UNIT_WARNING_STRING.format('Charges', 'elementary charge'))
            charge *= u.elementary_charge
        else:
            ensure_valid_dimensions(charge, u.elementary_charge)

        return charge

    @validator('mass_')
    def is_valid_mass(cls, mass):
        if mass is None:
            return None
        default_mass_units = u.gram /u.mol
        if not isinstance(mass, u.unyt_array):
            warnings.warn(UNIT_WARNING_STRING.format('Masses', 'g/mol'))
            mass *= default_mass_units
        else:
            ensure_valid_dimensions(mass, default_mass_units)
        return mass

    class Config:
        extra = 'forbid'

        fields = {
            'charge_': 'charge',
            'mass_': 'mass',
            'element_': 'element',
            'atom_type_': 'atom_type'
        }

        alias_to_fields = {
            'charge': 'charge_',
            'mass': 'mass_',
            'element': 'element_',
            'atom_type': 'atom_type_'
        }

        validate_assignment = True
