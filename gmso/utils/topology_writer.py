import warnings

from lxml import etree

from sympy import sympify
import unyt as u

from gmso.utils.ff_utils import validate

ATOM_TYPE_ATTRIBS = ['name',
                     'charge',
                     'mass',
                     'definition',
                     'description',
                     'atomclass',
                     'doi',
                     'overrides',
                     'definition',
                     'description']
BOND_TYPE_ATTRIBS = ['name', 'type1', 'type2']
ANGLE_TYPE_ATTRIBS = ['name', 'type1', 'type2', 'type3']
DIHEDRAL_TYPE_ATTRIBS = ['name', 'type1', 'type2', 'type3', 'type4']


def _parse_parameters_units(collection, check_consistency=False):
    for expression in collection:
        expression_symbols = sympify(expression).free_symbols
        param_units = {}
        for a_type in collection[expression]:
            parameters = expression_symbols - a_type.independent_variables
            for param in parameters:
                param_units[str(param)] = str(a_type.parameters[str(param)].units)
            if not check_consistency:
                break
    return param_units


def _populate_attribs(etree_el, obj, attribs_list):
    if hasattr(obj, 'member_types'):
        for i, member in enumerate(getattr(obj, 'member_types')):
            setattr(obj, 'type{}'.format(i+1), obj.member_types[i])

    for _attrib in attribs_list:
        val = getattr(obj, _attrib)
        if type(val) == u.unyt_quantity:
            val = str(val.value)
        etree_el.attrib[_attrib] = str(val)

    if len(obj.parameters) > 0:
        parameters = etree.SubElement(etree_el, 'Parameters')
        for name, val in obj.parameters.items():
            parameter = etree.SubElement(parameters, 'Parameter')
            parameter.attrib['name'] = name
            if type(val) == u.unyt_quantity:
                parameter.attrib['value'] = str(val.value)
            else:
                parameter.attrib['value'] = str(val)
    if hasattr(obj, 'member_types'):
        for i, member in enumerate(getattr(obj, 'member_types')):
            delattr(obj, 'type{}'.format(i+1))


def _assign_expressions(parent_tag, tag_name, expression):
    this_tag = etree.SubElement(parent_tag, tag_name)
    this_tag.attrib['expression'] = expression
    return this_tag


def _assign_parameters(params_unit_dict, parent_tag):
    for param, unit in params_unit_dict.items():
        param_unit_el = etree.SubElement(parent_tag, 'ParametersUnitDef')
        param_unit_el.attrib['parameter'] = param
        param_unit_el.attrib['unit'] = unit


def _create_element_tree(types_dict, parent_el, types):
    _type = types[0:-1]
    attributes = {'AtomType': ATOM_TYPE_ATTRIBS,
                  'BondType': BOND_TYPE_ATTRIBS,
                  'AngleType': ANGLE_TYPE_ATTRIBS,
                  'DihedralType': DIHEDRAL_TYPE_ATTRIBS}
    for expr in types_dict:
        # Create AtomTypes
        types_grp = _assign_expressions(parent_el, types, expr)
        params_unit_dict = _parse_parameters_units(types_dict)
        _assign_parameters(params_unit_dict, types_grp)
        for core_type in types_dict[expr]:
            type_el = etree.SubElement(types_grp, _type)
            _populate_attribs(type_el, core_type, attributes[_type])


class TopologyXMLWriter:
    @staticmethod
    def write(**kwargs):
        """Write a topology xml based on the key-word arguments"""
        # Create the root element
        ff_el = etree.Element('ForceField')
        ff_el.attrib['name'] = kwargs.get('name', 'ForceField')
        ff_el.attrib['version'] = kwargs.get('version', '1.0.0')

        # Create FFMetaData
        ff_meta = etree.SubElement(ff_el, 'FFMetaData')
        # Add Scaling Factors
        scaling_factors = kwargs.get('scaling_factors', {})
        ff_meta.attrib['electrostatics14Scale'] = str(scaling_factors.get('electrostatics14Scale', 1.0))
        ff_meta.attrib['nonBonded14Scale'] = str(scaling_factors.get('nonBonded14Scale', 1.0))

        # Create Units
        ff_units = etree.SubElement(ff_meta, 'Units')
        for unit_name, unit_val in kwargs.get('units', {}).items():
            ff_units.attrib[unit_name] = repr(unit_val)

        # Check any atom_types
        atom_types_dict = kwargs.get('atom_types', {})
        bond_types_dict = kwargs.get('bond_types', {})
        angle_types_dict = kwargs.get('angle_types', {})
        dihedral_types_dict = kwargs.get('dihedral_types', {})

        _create_element_tree(atom_types_dict, ff_el, 'AtomTypes')
        _create_element_tree(bond_types_dict, ff_el, 'BondTypes')
        _create_element_tree(angle_types_dict, ff_el, 'AngleTypes')
        _create_element_tree(dihedral_types_dict, ff_el, 'DihedralTypes')

        ff_tree = etree.ElementTree(ff_el)
        validate(ff_tree)
        ff_tree.write(kwargs.get('filename', 'forcefield.xml'),
                      pretty_print=True,
                      xml_declaration=True,
                      encoding='utf-8')





