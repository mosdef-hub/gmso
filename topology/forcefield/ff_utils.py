import os

import unyt as u
from sympy import sympify
from lxml import etree

from topology.core.atom_type import AtomType
from topology.core.bond_type import BondType
from topology.core.angle_type import AngleType
from topology.core.dihedral_type import DihedralType
from topology.exceptions import ForceFieldParseError

__all__ = ['validate', 'parse_ff_metadata', 'parse_ff_atomtypes', 'parse_ff_connection_types']


def _parse_param_units(parent_tag):
    param_unit_dict = {}
    params_iter = parent_tag.getiterator('ParametersUnitDef')
    for param_unit in params_iter:
        param_unit_dict[param_unit.attrib['parameter']] = param_unit.attrib['unit']
    return param_unit_dict


def _parse_params_values(parent_tag, units_dict):
    # Tag of type Parameters can exist atmost once
    params_dict = {}
    for param in parent_tag.find('Parameters').getiterator('Parameter'):
        if param.attrib['name'] not in units_dict:
            raise ForceFieldParseError('Parameters {} with Unknown units found'.format(param.attrib['name']))
        param_name = param.attrib['name']
        param_unit = units_dict[param_name]
        param_value = u.unyt_quantity(float(param.attrib['value']), param_unit)
        params_dict[param_name] = param_value
    return params_dict


def _check_valid_atomtype_names(tag, ref_dict):
    at1 = tag.attrib.get('type1', None)
    at2 = tag.attrib.get('type2', None)
    at3 = tag.attrib.get('type3', None)
    at4 = tag.attrib.get('type4', None)

    member_types = list(filter(lambda x: x is not None, [at1, at2, at3, at4]))

    for member in member_types:
        if member not in ref_dict:
            raise ForceFieldParseError('AtomTypes {} not present in AtomTypes reference in the xml'.format(member))

    return member_types


def _parse_units(unit_tag):
    if unit_tag is None:
        unit_tag = {}
    units_map = {
        'energy': u.kcal / u.mol,
        'distance': u.nm,
        'mass': u.amu,
        'charge': u.coulomb
    }
    for attrib, val in unit_tag.items():
        units_map[attrib] = u.Unit(val)
    return units_map


def validate(xml_path, schema=None):
    """Validate a given xml file with a refrence schema"""
    if schema is None:
        schema_path = os.path.join(os.path.split(os.path.abspath(__file__))[0], 'schema', 'ff-topology.xsd')

    xml_doc = etree.parse(schema_path)
    xmlschema = etree.XMLSchema(xml_doc)
    ff_xml = etree.parse(xml_path)
    xmlschema.assertValid(ff_xml)


def parse_ff_metadata(element):
    metatypes = ['Units']
    parsers = {
        'Units': _parse_units
    }
    ff_meta = {}
    for metatype in element:
        if metatype.tag in metatypes:
            ff_meta[metatype.tag] = parsers[metatype.tag](metatype)
    return ff_meta


def parse_ff_atomtypes(atomtypes_el, ff_meta):
    """Given an xml element tree rooted at AtomType, traverse the tree to form a proper topology.core.AtomType"""
    atomtypes_dict = {}
    units_dict = ff_meta['Units']
    atom_types_expression = atomtypes_el.attrib['expression']
    param_unit_dict = _parse_param_units(atomtypes_el)

    # Parse all the atomTypes and create a new AtomType
    for atom_type in atomtypes_el.getiterator('AtomType'):
        ctor_kwargs = {
            'name': 'AtomType',
            'mass': 0.0 * u.elementary_charge,
            'expression': '4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
            'parameters': None,
            'independent_variables': None,
            'atomclass': '',
            'doi': '',
            'overrides': '',
            'definition': '',
            'description': '',
            'topology': None
        }

        if atom_types_expression:
            ctor_kwargs['expression'] = atom_types_expression

        for kwarg in ctor_kwargs.keys():
            ctor_kwargs[kwarg] = atom_type.attrib.get(kwarg, ctor_kwargs[kwarg])
        if isinstance(ctor_kwargs['mass'], str):
            ctor_kwargs['mass'] = u.unyt_quantity(float(ctor_kwargs['mass']), units_dict['mass'])
        if isinstance(ctor_kwargs['overrides'], str):
            ctor_kwargs['overrides'] = set(ctor_kwargs['overrides'].split(','))
        params_dict = _parse_params_values(atom_type, param_unit_dict)
        if not ctor_kwargs['parameters']:
            ctor_kwargs['parameters'] = params_dict
        valued_param_vars = set(sympify(param) for param in params_dict.keys())
        ctor_kwargs['independent_variables'] = sympify(atom_types_expression).free_symbols - valued_param_vars
        this_atom_type = AtomType(**ctor_kwargs)
        atomtypes_dict[this_atom_type.name] = this_atom_type
    return atomtypes_dict


TAG_TO_CLASS_MAP = {
    'BondType': BondType,
    'AngleType': AngleType,
    'DihedralType': DihedralType
}


def parse_ff_connection_types(connectiontypes_el, atomtypes_dict, child_tag='BondType'):
    """Given an XML etree Element rooted at BondTypes, parse the XML to create topology.core.AtomTypes,"""
    connectiontypes_dict = {}
    connectiontype_expression = connectiontypes_el.attrib['expression']
    param_unit_dict = _parse_param_units(connectiontypes_el)

    # Parse all the bondTypes and create a new BondType
    for connection_type in connectiontypes_el.getiterator(child_tag):
        ctor_kwargs = {
            'name': child_tag,
            'expression': '0.5 * k * (r-r_eq)**2',
            'parameters': None,
            'independent_variables': None,
            'member_types': None
        }
        if connectiontype_expression:
            ctor_kwargs['expression'] = connectiontype_expression

        for kwarg in ctor_kwargs.keys():
            ctor_kwargs[kwarg] = connection_type.attrib.get(kwarg, ctor_kwargs[kwarg])

        ctor_kwargs['member_types'] = _check_valid_atomtype_names(connection_type, atomtypes_dict)
        if not ctor_kwargs['parameters']:
            ctor_kwargs['parameters'] = _parse_params_values(connection_type, param_unit_dict)
        valued_param_vars = set(sympify(param) for param in ctor_kwargs['parameters'].keys())
        ctor_kwargs['independent_variables'] = sympify(connectiontype_expression).free_symbols - valued_param_vars

        this_conn_type = TAG_TO_CLASS_MAP[child_tag](**ctor_kwargs)
        connectiontypes_dict[this_conn_type.name] = this_conn_type
    return connectiontypes_dict
