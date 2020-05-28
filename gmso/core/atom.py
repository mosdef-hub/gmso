import warnings
from typing import Optional, Union

import unyt as u

from pydantic import validator

from gmso.core.atom_type import AtomType
from gmso.core.element import Element
from gmso.abc.abstract_site import Site


class Atom(Site):
    """An interaction site object in the topology hierarchy.

    Site is the object that represents any general interaction site in a molecular simulation.
    Sites have been designed to be as general as possible, making no assumptions about representing atoms or beads, or having mass or charge.
    That is, a Site can represent an atom in an atomistic system, a bead in a coarse-grained system, and much more.

    Parameters
    ----------
    name : str, optional, default='Site'
       Name of the site
    position : unyt array or numpy array or list, optional, default=None
       The position of the site in Cartesian space.
       If a unyt array is not passed, units are assumed to be in 'nm'.
    charge : unyt quantity or float, optional, default=None
       The charge of the site.
       Unyt quantities are converted to units of elementary charge, float values are assumed to be in units of elementary charge.  
       If no value is passed, site attempts to grab a charge from site.atom_type.
    mass : unyt quantity or float, optional, default=None
       The mass of the site.  
       Unyt quantities are converted to units of g/mol, float values are assumed to be in units of g/mol.  
       If no value is passed, site attempts to grab a mass from site.atom_type.
    element : 'Element' object, optional, default=None
       The element of the site represented by the `Element` object.  
       See `element.py` for more information.
    atom_type : 'AtomType' object, optional, default=None
       The atom type of the site containing functional forms, interaction parameters, and other properties such as mass and charge.  
       See `atom_type.py` for more information.

    """

    charge_: Optional[Union[u.unyt_quantity, float]] = None
    mass_: Optional[Union[u.unyt_quantity, float]] = None
    element_: Optional[Element] = None
    atom_type_: Optional[AtomType] = None

    @property
    def charge(self) -> Union[u.unyt_quantity, None]:
        if self.charge_:
            return self.charge_
        elif self.atom_type_:
            return self.atom_type_.charge
        return None

    @property
    def mass(self) -> Union[u.unyt_quantity, None]:
        if self.mass_:
            return self.mass_
        elif self.atom_type:
            return self.atom_type.mass
        return None

    @property
    def element(self):
        return self.element_

    @property
    def atom_type(self) -> Union[AtomType, None]:
        return self.__dict__['atom_type_']

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
