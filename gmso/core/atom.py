import warnings
from typing import Optional, Union

import unyt as u

from pydantic import Field, validator

from gmso.core.atom_type import AtomType
from gmso.core.element import Element
from gmso.abc.abstract_site import Site


class Atom(Site):
    __base_doc__ = """An atom represents a single element association in a topology.
    Atom is the object that represents any general interaction Atom in a molecular simulation.
    
    Notes
    -----
    Atoms have all the attributes inherited from the base Site class
    
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

    @validator('charge_')
    def is_valid_charge(cls, charge):
        if not isinstance(charge, u.unyt_array):
            warnings.warn("Charges are assumed to be elementary charge")
            charge *= u.elementary_charge

        elif charge.units.dimensions != u.elementary_charge.units.dimensions:
            warnings.warn("Charges are assumed to be elementary charge")
            charge = charge.value * u.elementary_charge

        return charge

    @validator('element_')
    def is_valid_element(cls, element):
        if element is not None:
            assert isinstance(element, Element), 'Attribute element must be of type element'
        return element

    @validator('mass_')
    def is_valid_mass(cls, mass):
        if not isinstance(mass, u.unyt_array):
            warnings.warn("Masses are assumed to be g/mol")
            mass *= u.gram / u.mol

        elif mass.units.dimensions != (u.gram / u.mol).units.dimensions:
            warnings.warn("Masses are assumed to be g/mol")
            mass = mass.value * u.gram / u.mol
        return mass

    class Config(Site.Config):
        extra = 'forbid'
        fields = {
            'charge_': 'charge',
            'mass_': 'mass',
            'element_': 'element',
            'atom_type_': 'atom_type'
        }
        # For __setattr__
        Site.Config.alias_to_fields.update({
            'charge': 'charge_',
            'mass': 'mass_',
            'element': 'element_',
            'atom_type': 'atom_type_'
        })
        validate_assignment = True
