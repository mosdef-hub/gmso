from collections import ChainMap

from lxml import etree

from gmso.utils.ff_utils import (validate,
                                 parse_ff_metadata,
                                 parse_ff_atomtypes,
                                 parse_ff_connection_types)


class ForceField(object):
    """A generic implementation of the forcefield class.

    The ForceField class is one of the core data structures in gmso, which is
    used to hold a collection of gmso.core.Potential subclass objects along with some
    metadata to represent a forcefield. The forcefield object can be applied
    to any gmso.Topology which has effects on its Sites, Bonds, Angles and Dihedrals.

    Parameters
    ----------
    name : str
        Name of the forcefield, default 'ForceField'
    version : str
        a cannonical semantic version of the forcefield, default 1.0.0

    Attributes
    ----------
    name : str
        Name of the forcefield
    version : str
        Version of the forcefield
    atom_types : dict
        A collection of atom types in the forcefield
    bond_types : dict
        A collection of bond types in the forcefield
    angle_types : dict
        A collection of angle types in the forcefield
    dihedral_types : dict
        A collection of dihedral types in the forcefield
    units : dict
        A collection of unyt.Unit objects used in the forcefield
    scaling_factors : dict
        A collection of scaling factors used in the forcefield

    See Also
    --------
    gmso.ForceField.from_xml : A class method to create forcefield object from XML files

    """
    def __init__(self, xml_loc=None):
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
                atomclass_group = atomclass_dict.get(atom_type.atomclass, [])
                atomclass_group.append(atom_type)
                atomclass_dict[atom_type.atomclass] = atomclass_group
        return atomclass_dict

    @classmethod
    def from_xml(cls, xml_locs):
        """Create a gmso.Forcefield object from XML File(s)

        This class method creates a ForceFiled object from the reference
        XML file. This method takes in a single or collection of XML files
        with information about gmso.AtomTypes, gmso.BondTypes, gmso.AngleTypes
        and gmso.DihedralTypes to create the ForceField object.

        Parameters
        ----------
        xml_locs : str or iterable of str
            string or iterable of strings containing the forcefield XML locations

        Returns
        --------
        forcefield : gmso.ForceField
            A gmso.Forcefield object with a collection of Potential objects
            created using the information in the XML file
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

        atom_types_dict = ChainMap()
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
        ff.atom_types = atom_types_dict.maps[0]
        ff.bond_types = bond_types_dict
        ff.angle_types = angle_types_dict
        ff.dihedral_types = dihedral_types_dict
        return ff
