from enum import IntEnum

from lxml import etree

from topology.forcefield.ff_utils import validate

__all__ = ['XMLS', 'XMLConverter']


class XMLS(IntEnum):
    TOPOLOGY = 1
    OPENMM = 2


class TopologyXMLWriter(object):

    @staticmethod
    def write(consolidated_dict, dst):
        """Write a consolidated into a XML File"""
        ff_dict = consolidated_dict.get('forcefield')
        # Write The root
        ff_el = etree.Element('ForceField')
        ff_el.attrib['name'] = ff_dict.get('name', 'ForceField')
        ff_el.attrib['version'] = ff_dict.get('version', '1.0.0')

        # Write MetaData and units
        ff_meta_el = etree.SubElement(ff_el, 'FFMetaData')
        ff_scaling_factors = ff_dict.get('scaling_factors', {
            'electrostatics14Scale': 1.0,
            'nonBonded14Scale': 1.0
        })
        for key, value in ff_scaling_factors.items():
            ff_meta_el.attrib[key] = str(value)

        ff_units = ff_dict.get('units', None)
        if ff_units:
            units = etree.SubElement(ff_meta_el, 'Units')
            for key, val in ff_units.items():
                units.attrib[key] = str(val)

        atom_types_list = consolidated_dict.get('atom_types')
        for atom_types in atom_types_list:
            atom_types_el = etree.SubElement(ff_el, 'AtomTypes')
            atom_types_el.attrib['expression'] = atom_types.get('expression', '')
            for parameters in atom_types.get('parameters_units_def', []):
                parameters_units_def_el = etree.SubElement(atom_types_el, 'ParametersUnitDef')
                for param, unit in parameters.items():
                    parameters_units_def_el.attrib['parameter'] = param
                    parameters_units_def_el.attrib['unit'] = unit
            for atom_type in atom_types.get('collection', []):
                atom_type_el = etree.SubElement(atom_types_el, 'AtomType')
                for key, val in atom_type.items():
                    if key != 'parameters':
                        atom_type_el.attrib[key] = val
                    else:
                        parameters_el = etree.SubElement(atom_type_el, 'Parameters')
                        for parameter in atom_type[key]:
                            parameter_el = etree.SubElement(parameters_el, 'Parameter')
                            parameter_el.attrib['name'] = parameter['name']
                            parameter_el.attrib['value'] = str(parameter['value'])

        ff_tree = etree.ElementTree(ff_el)
        validate(ff_tree)
        ff_tree.write(dst, pretty_print=True, xml_declaration=True, encoding='utf-8')


class _XMLParser(object):
    def __init__(self, xml_file):
        self.ff_tree = etree.parse(xml_file)
        self.topo_writer = TopologyXMLWriter

    def parse_to_topo_dict(self):
        pass

    def write_topo_xml(self, consolidated_dict, dst):
        self.topo_writer.write(consolidated_dict, dst)


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
            'dihedral_types': self._parse_dihedral_types()
        }

    def _parse_forcefield(self):
        return {
            'name': self.ff_tree.getroot().tag,
            'version': '1.0.0',
            'units': {
                'distance': 'nm',
                'time': 'ps',
                'mass': 'amu',
                'charge': 'qp',
                'angle': 'rad',
                'energy': 'KJ/mol',
                'temperature': 'K'
            },
            'scaling_factors': self._parse_scaling_factors()
        }

    def _parse_scaling_factors(self):
        non_bonded_forces_el = self.ff_tree.findall('NonbondedForce')
        scaling_factors = {
            'electrostatics14Scale': 1.0,
            'nonBonded14Scale': 1.0
        }
        for el in non_bonded_forces_el:
            scaling_factors['electrostatics14Scale'] = float(el.get('coulomb14scale', 1.0))
            scaling_factors['nonBonded14Scale'] = float(el.get('lj14scale', 1.0))

        return scaling_factors

    def _parse_atom_types(self):
        atom_types_el = self.ff_tree.findall('AtomTypes')
        non_bonded_forces_el = self.ff_tree.findall('NonbondedForce')
        atom_types_list = []
        for atom_types in atom_types_el:
            for atom_type in atom_types.getiterator('Type'):
                # Get the atom_type properties here
                this_atom_type = {
                    'name': atom_type.get('name', 'AtomType'),
                    'atomclass': atom_type.get('class', ''),
                    'element': atom_type.get('element', ''),
                    'mass': atom_type.get('mass', '1.0'),
                    'charge': atom_type.get('charge', '1.0'),
                    'definition': atom_type.get('def', '')
                }
                # Check to see if this atom_type is somewhere in the
                # <NonbondedForce>
                for el in non_bonded_forces_el:
                    atom_type_tag = el.find('.//Atom[@type="{}"]'.format(this_atom_type['name']))
                    this_atom_type['parameters'] = [{'name': 'sigma', 'value': float(atom_type_tag.get('sigma', 1.0))},
                                                    {'name': 'e0', 'value': 8.8542e-12},
                                                    {'name': 'q', 'value': float(atom_type_tag.get('charge', 1.0))},
                                                    {'name': 'ep', 'value': float(atom_type_tag.get('epsilon', 1.0))}]
                    this_atom_type['charge'] = atom_type_tag.get('charge', 1.0)

                atom_types_list.append(this_atom_type)
        return [{
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
                    'unit': 'qp'
                }
            ],
            'collection': atom_types_list
        }]

    def _parse_bond_types(self):
        # harmonic_bond_force_el = self.ff_tree.findall('HarmonicBondForce')
        # for bond_types in harmonic_bond_force_el:
        pass

    def _parse_angle_types(self):
        pass

    def _parse_dihedral_types(self):
        pass


class XMLConverter(object):
    """A generic base class to convert external XML files to
    topology XML schema and vice versa. This class provides a
    a single static method that handles all the cases
    of external xml conversion.
    """
    CONVERTERS = {
        XMLS.OPENMM: _OpenMMXMLParser
    }

    @staticmethod
    def convert(src, dst, src_format=XMLS.OPENMM):
        """Convert a source forcefield XML(from any external standard to )"""
        converter = XMLConverter.CONVERTERS[src_format](src)
        topo_dict = converter.parse_to_topo_dict()
        converter.write_topo_xml(topo_dict, dst)


