from enum import IntEnum

from lxml import etree

__all__ = ['XMLS', 'XMLConverter']


class XMLS(IntEnum):
    TOPOLOGY = 1
    OPENMM = 2


class _XMLParser(object):
    def __init__(self, xml_file):
        self.ff_tree = etree.parse(xml_file)

    def parse_to_topo_dict(self):
        pass


class _OpenMMXMLParser(_XMLParser):
    def __init__(self, xml_file):
        super(_OpenMMXMLParser, self).__init__(xml_file)

    def parse_to_topo_dict(self):
        """Parse the contents of the OpenMM Forcefield XML into a dictionary as
            {
                forcefield: object representing forcefield name/version(optional)
                atom_types: All possible mappings from openMM XML tags to AtomTypes
                bond_types: All possible mappings from openMM XML tags to BondTypes
                angle_types: All possible mappings from openMM XML tags to AngleTypes
                dihedral_types: All possible mappings from openMM XML tags to DihedralTypes
            }
        """
        return {
            'forcefield': self._parse_forcefield(),
            'atom_types': self._parse_atom_types(),
            'bond_types': self._parse_bond_types(),
            'angle_types': self._parse_angle_types(),
        }

    def _parse_forcefield(self):
        return {
            'name': self.ff_tree.root.tag,
            'version': '1.0.0',
            'ff_meta': {
                'units': {
                    'distance': 'nm',
                    'time': 'ps',
                    'mass': 'amu',
                    'charge': 'qp',
                    'angle': 'rad',
                    'energy': 'KJ/mol'
                }
            }
        }

    def _parse_atom_types(self):
        atom_types_el = self.ff_tree.findall('AtomTypes')
        non_bonded_forces_el = self.ff_tree.findall('NonbondedForce')
        atom_types_list = []
        for atom_types in atom_types_el:
            for atom_type in atom_types.getiterator('Type'):
                this_atom_type = {
                    'name': atom_type.get('name', 'AtomType'),
                    'atomclass': atom_type.get('class', ''),
                    'element': atom_type.get('element', ''),
                    'mass': atom_type.get('mass', '1.0'),
                    'charge': atom_type.get('charge', '1.0'),
                    'definition': atom_type.get('def', '')
                }

                atom_types_list.append(this_atom_type)
        return {
            'expression': 'ep * ((sigma/r)**12 - (sigma/r)**6) + q / (e0 * r)',
            'parameters_units_def': [
                {
                    'parameter': 'ep',
                    'unit': 'KJ/mol',
                },
                {
                    'parameter': 'sigma',
                    'unit': 'nm'
                },
                {
                    'parameter': 'e0',
                    'unit': 'A**2*s**4/(kg*m**3)'
                },
                {
                    'parameter': 'q',
                    'unit': 'coulomb'
                }
            ],
            'collection': atom_types_list
        }

    def _parse_bond_types(self):
        pass

    def _parse_angle_types(self):
        pass


class XMLConverter(object):
    """A generic base class to convert external XML files to
    topology XML schema and vice versa. This class provides a
    a single static method that handles all the cases
    of external xml conversion.
    """
    CONVERTERS = {
        XMLS.TOPOLOGY: _XMLParser,
        XMLS.OPENMM: _OpenMMXMLParser
    }
    @staticmethod
    def convert(src, dst, src_format=XMLS.OPENMM, dst_format=XMLS.OPENMM):
        """Convert a source forcefield XML(from any external standard to )"""
        converter = XMLConverter.CONVERTERS[src_format](src)
        converter.parse_to_topo_dict()





