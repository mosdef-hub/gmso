import xml.etree.ElementTree as ET

from topology.core.atom_type import AtomType
from .ff_utils import mass_to_unyt
class ForceField(object):
    """A Generic implementation of the forcefield class
    A forcefield class contains different collection of
    core type members.
    Parameters:
    ----------
    """
    def __init__(self, name='ForceField'):
        """Initialize a new ForceField"""
        self.name = name
        self.atom_types = {}
        self.bond_types = {}
        self.angle_types = {}
        self.dihedral_types = {}

    def __repr__(self):
        descr = list('<Forcefield ')
        descr.append(self.name + ' ')
        descr.append('{:d} AtomTypes, '.format(len(self.atom_types)))
        descr.append('{:d} BondTypes, '.format(len(self.bond_types)))
        descr.append('{:d} AngleTypes, '.format(len(self.angle_types)))
        descr.append('id: {}>'.format(id(self)))
        return ''.join(descr)

    @classmethod
    def from_xml(cls, xml_locs):
        """Create a forcefield object from a XML File
        Parameters:
        -----------
        xml_locs: (str) or iter(str), string or iteratble of strings
                  containing the forcefield XML locations

        Returns:
        --------
        topology.forcefield.ForceField object, containing all the information
            from the ForceField File
        """
        atom_type_map = {}
        atom_types = {}
        if not hasattr(xml_locs, '__iter__'):
            xml_locs = [].append(xml_locs)
        if isinstance(xml_locs, str):
            xml_locs = [xml_locs]
        for loc in xml_locs:
            ff_tree = ET.parse(loc)
            for atom_types in ff_tree.findall('AtomTypes'):
                for atom_type in atom_types:
                    atom_type_props = atom_type.attrib
                    this_atom_type = AtomType(
                        name=atom_type_props['name'],
                        mass=mass_to_unyt(atom_type_props['mass'])
                    )
                    if this_atom_type in atom_type_map:
                        atom_types[atom_type_props['name']] = atom_type_map[this_atom_type]
                    else:
                        atom_type_map[this_atom_type] = this_atom_type
