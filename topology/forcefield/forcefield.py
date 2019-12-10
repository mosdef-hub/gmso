from lxml import etree

from topology.forcefield.ff_utils import validate, parse_ff_metadata, parse_ff_atomtype


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

        if not hasattr(xml_locs, '__iter__'):
            xml_locs = [].append(xml_locs)
        if isinstance(xml_locs, str):
            xml_locs = [xml_locs]
        for loc in xml_locs:
            validate(loc)
            ff_tree = etree.parse(loc)
            ff_meta_tree = ff_tree.find('FFMetaData')
            ff_meta_map = parse_ff_metadata(ff_meta_tree)
            ff_atomtypes_list = ff_tree.findall('AtomTypes')
            atom_types_dict = {}
            for atom_types in ff_atomtypes_list:
                atom_types_dict.update(parse_ff_atomtype(atom_types, ff_meta_map))

        ff = cls()
        ff.atom_types = atom_types_dict
        return ff
