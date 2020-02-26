from lxml import etree

from gmso.forcefield.ff_utils import (validate,
                                          parse_ff_metadata,
                                          parse_ff_atomtypes,
                                          parse_ff_connection_types)


class ForceField(object):
    """A generic implementation of the forcefield class.

    A forcefield class contains different collection of
    core type members.

    Parameters
    ----------
    name: (str), Name of the forcefield, default 'ForceField'
    version: (str), a cannonical semantic version of the forcefield, default 1.0.0

    Attributes
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

        Parameters
        ----------
        name: str, name of the forcefield, default 'ForceField', usage optional
        version: str, version of the forcefield, default 'ForceField', usage optional
        xml_locs: iterable, list of GMSO's Forcefield XML Files, default None, usage optional
        """
        if xml_loc is not None:
            ff = ForceField.from_xml(xml_loc)
            self.name = ff.name
            self.version = ff.version
            self.atom_types = ff.atom_types
            self.bond_types = ff.bond_types
            self.angle_types = ff.angle_types
            self.dihedral_types = ff.dihedral_types
            self.scaling_factors = ff.scaling_factors
            self.units = ff.units
        else:
            self.name = 'ForceField'
            self.version = '1.0.0'
            self.atom_types = {}
            self.bond_types = {}
            self.angle_types = {}
            self.dihedral_types = {}
            self.scaling_factors = {}
            self.units = {}

    def __repr__(self):
        descr = list('<Forcefield ')
        descr.append(self.name + ' ')
        descr.append('{:d} AtomTypes, '.format(len(self.atom_types)))
        descr.append('{:d} BondTypes, '.format(len(self.bond_types)))
        descr.append('{:d} AngleTypes, '.format(len(self.angle_types)))
        descr.append('{:d} DihedralTypes, '.format(len(self.dihedral_types)))
        descr.append('id: {}>'.format(id(self)))
        return ''.join(descr)

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
        gmso.forcefield.ForceField object, containing all the information
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
        ff.scaling_factors = ff_meta_map['scaling_factors']
        ff.units = ff_meta_map['Units']
        ff.atom_types = atom_types_dict
        ff.bond_types = bond_types_dict
        ff.angle_types = angle_types_dict
        ff.dihedral_types = dihedral_types_dict
        return ff
