import re

from lxml import etree

from topology.forcefield.ff_utils import (validate,
                                          parse_ff_metadata,
                                          parse_ff_atomtypes,
                                          parse_ff_connection_types,
                                          DICT_KEY_SEPARATOR)
from topology.exceptions import ForceFieldError


class ForceField(object):
    """A Generic implementation of the forcefield class
    A forcefield class contains different collection of
    core type members.
    Parameters:
    ----------
    name: (str), Name of the forcefield, default 'ForceField'
    version: (str), a cannonical semantic version of the forcefield, default 1.0.0

    Attributes:
    -----------
    name: (str), Name of the forcefield
    version: (str), Version of the forcefield
    atom_types: (dict), A collection of atom types in the forcefield
    bond_types: (dict), A collection of bond types in the forcefield
    angle_types: (dict), A collection of angle types in the forcefield
    dihedral_types: (dict), A collection of dihedral types in the forcefield
    """

    def __init__(self, xml_loc=None):
        """Initialize a new ForceField
        Parameters:
        ------------
        name: str, name of the forcefield, default 'ForceField', usage optional
        version: str, version of the forcefield, default 'ForceField', usage optional
        xml_locs: iterable, list of topology's Forcefield XML Files, default None, usage optional
        """
        if xml_loc is not None:
            ff = ForceField.from_xml(xml_loc)
            self.name = ff.name
            self.version = ff.version
            self.atom_types = ff.atom_types
            self.bond_types = ff.bond_types
            self.angle_types = ff.angle_types
            self.dihedral_types = ff.dihedral_types
        else:
            self.name = 'ForceField'
            self.version = '1.0.0'
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
        descr.append('{:d} DihedralTypes, '.format(len(self.dihedral_types)))
        descr.append('id: {}>'.format(id(self)))
        return ''.join(descr)

    def __getitem__(self, item=''):
        """Get AtomType, BondType, AngleType or DihedralType in a ForceField with support for wildcards"""
        splitted_items = item.split(DICT_KEY_SEPARATOR)
        keys_map = {
            1: self.atom_types,
            2: self.bond_types,
            3: self.angle_types,
            4: self.dihedral_types
        }
        ordinals = ['1st', '2nd', '3rd', '4th']

        if len(splitted_items) > 4:
            raise ForceFieldError('Error: You provided a key of length {}. Keys longer than 4 items are not supported '
                                  'by the forcefield class.')

        # Pass 1: Direct Lookup
        try:
            return keys_map[len(splitted_items)][item]
        # Pass 3: WildCard Map
        except KeyError:
            class_names = ['AtomType', 'BondType', 'AngleType', 'DihedralType']
            wildcard_idxes = []
            wildcard = '*'

            for idx, item in enumerate(splitted_items):
                if item == wildcard:
                    if idx != 0 and idx != len(splitted_items) - 1:
                        print(idx, len(splitted_items))
                        raise ForceFieldError('Error. Wildcard for {} to replace {} element is not supported.'
                                              .format(class_names[len(splitted_items) - 1], ordinals[idx]))
                    wildcard_idxes.append(idx)

            if len(wildcard_idxes) == 0:
                # Try a reverse lookup
                try:
                    return keys_map[len(splitted_items)][DICT_KEY_SEPARATOR.join(reversed(splitted_items))]
                except KeyError:
                    raise KeyError('No Matching {0} for {1} found in the Forcefield.'
                                   .format(class_names[len(splitted_items) - 1],
                                           DICT_KEY_SEPARATOR.join(splitted_items)))
            possible_keys = keys_map[len(splitted_items)].keys()

            pattern_list = []
            for item in splitted_items:
                if item == wildcard:
                    pattern_list.append('.*')
                else:
                    pattern_list.append(item)
            for key in possible_keys:
                if re.match(DICT_KEY_SEPARATOR.join(pattern_list), key):
                    return self.__getitem__(key)
                if re.match(DICT_KEY_SEPARATOR.join(reversed(pattern_list)), key):
                    return self.__getitem__(key)
            raise KeyError("{0} {1} doesn't exist in this ForceField".format(class_names[len(splitted_items) - 1],
                                                                             DICT_KEY_SEPARATOR.join(splitted_items)))

    @property
    def atom_class_groups(self):
        """Return a dictionary of atomClasses in the Forcefield"""
        atom_types = self.atom_types.values()
        atomclass_dict = {}
        for atom_type in atom_types:
            if atom_type.atomclass is not None:
                this_atomtype_class_group = atomclass_dict.get(atom_type.atomclass, [])
                this_atomtype_class_group.append(atom_type)
        return atomclass_dict

    @classmethod
    def from_xml(cls, xml_locs):
        """Create a forcefield object from a XML File
        Parameters:
        -----------
        xml_locs: (str) or iter(str), string or iterable of strings
                  containing the forcefield XML locations

        Returns:
        --------
        topology.forcefield.ForceField object, containing all the information
            from the ForceField File
        """

        if not hasattr(xml_locs, '__iter__'):
            xml_locs = [].append(xml_locs)
        if isinstance(xml_locs, str):
            xml_locs = [xml_locs]

        versions = []
        names = []
        ff_atomtypes_list = []
        ff_bondtypes_list = []
        ff_angletypes_list = []
        ff_dihedraltypes_list = []

        atom_types_dict = {}
        bond_types_dict = {}
        angle_types_dict = {}
        dihedral_types_dict = {}

        for loc in xml_locs:
            validate(loc)
            ff_tree = etree.parse(loc)
            ff_el = ff_tree.getroot()
            versions.append(ff_el.attrib['version'])
            names.append(ff_el.attrib['name'])
            ff_meta_tree = ff_tree.find('FFMetaData')

            if ff_meta_tree is not None:
                ff_meta_map = parse_ff_metadata(ff_meta_tree)

            ff_atomtypes_list.extend(ff_tree.findall('AtomTypes'))
            ff_bondtypes_list.extend(ff_tree.findall('BondTypes'))
            ff_angletypes_list.extend(ff_tree.findall('AngleTypes'))
            ff_dihedraltypes_list.extend(ff_tree.findall('DihedralTypes'))

        # Consolidate AtomTypes
        for atom_types in ff_atomtypes_list:
            atom_types_dict.update(parse_ff_atomtypes(atom_types, ff_meta_map))

        # Consolidate BondTypes
        for bond_types in ff_bondtypes_list:
            bond_types_dict.update(parse_ff_connection_types(bond_types,
                                                             atom_types_dict,
                                                             child_tag='BondType'))

        # Consolidate AngleTypes
        for angle_types in ff_angletypes_list:
            angle_types_dict.update(parse_ff_connection_types(angle_types,
                                                              atom_types_dict,
                                                              child_tag='AngleType'))

        # Consolidate DihedralTypes
        for dihedral_types in ff_dihedraltypes_list:
            dihedral_types_dict.update(parse_ff_connection_types(dihedral_types,
                                                                 atom_types_dict,
                                                                 child_tag='DihedralType'))

        ff = cls()
        ff.name = names[0]
        ff.version = versions[0]
        ff.atom_types = atom_types_dict
        ff.bond_types = bond_types_dict
        ff.angle_types = angle_types_dict
        ff.dihedral_types = dihedral_types_dict
        return ff
