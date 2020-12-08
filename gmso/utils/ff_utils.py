import os
import re
from collections import ChainMap

import unyt as u
from sympy import sympify
from lxml import etree

from gmso.core.atom_type import AtomType
from gmso.core.bond_type import BondType
from gmso.core.angle_type import AngleType
from gmso.core.dihedral_type import DihedralType
from gmso.core.improper_type import ImproperType
from gmso.exceptions import ForceFieldParseError, ForceFieldError

__all__ = ['validate',
           'parse_ff_metadata',
           'parse_ff_atomtypes',
           'parse_ff_connection_types',
           'DICT_KEY_SEPARATOR']

DICT_KEY_SEPARATOR = '~'

# Create a dictionary of units
_unyt_dictionary = {}
for name, item in vars(u).items():
    if isinstance(item, u.Unit) or isinstance(item, u.unyt_quantity):
        _unyt_dictionary.update({name: item})


def _check_valid_string(type_str):
    if DICT_KEY_SEPARATOR in type_str:
        raise ForceFieldError('Please do not use {} in type string'.format(DICT_KEY_SEPARATOR))


def _parse_param_units(parent_tag):
    param_unit_dict = {}
    params_iter = parent_tag.getiterator('ParametersUnitDef')
    for param_unit in params_iter:
        param_unit_dict[param_unit.attrib['parameter']] = _parse_unit_string(param_unit.attrib['unit'])
    return param_unit_dict


def _parse_params_values(parent_tag, units_dict, child_tag, expression=None):
    # Tag of type Parameters can exist atmost once
    params_dict = {}
    if parent_tag.find('Parameters') is None:
        return params_dict
    for param in parent_tag.find('Parameters').getiterator('Parameter'):
        if param.attrib['name'] not in units_dict:
            raise ForceFieldParseError('Parameters {} with Unknown units found'.format(param.attrib['name']))
        param_name = param.attrib['name']
        param_unit = units_dict[param_name]
        param_value = u.unyt_quantity(float(param.attrib['value']), param_unit)
        params_dict[param_name] = param_value
    param_ref_dict = units_dict
    if child_tag == 'DihedralType':
        if not expression:
            raise ForceFieldError('Cannot consolidate parameters without an expression')
        _consolidate_params(params_dict, expression)
        param_ref_dict = _consolidate_params(units_dict, expression, update_orig=False)

    for param in param_ref_dict:
        if param not in params_dict:
            raise ForceFieldParseError(
                'Parameter {} is in units but cannot be found in parameters list'.format(param))
    return params_dict


def _consolidate_params(params_dict, expression, update_orig=True):
    to_del = []
    new_dict = {}
    match_string = '|'.join(str(symbol) for symbol in sympify(expression).free_symbols)
    for param in params_dict:
        match = re.match(r"({0})([0-9]+)".format(match_string), param)
        if match:
            new_dict[match.groups()[0]] = new_dict.get(match.groups()[0], [])
            new_dict[match.groups()[0]].append(params_dict[param])
            to_del.append(param)
    if update_orig:
        for key in to_del:
            del params_dict[key]
        params_dict.update(new_dict)
    else:
        return new_dict


def _check_atomtype_existence(tag, ref_dict):
    at1 = tag.attrib.get('type1', tag.attrib.get('class1', None))
    at2 = tag.attrib.get('type2', tag.attrib.get('class2', None))
    at3 = tag.attrib.get('type3', tag.attrib.get('class3', None))
    at4 = tag.attrib.get('type4', tag.attrib.get('class4', None))

    atom_classes_set = set(value.atomclass for value in ref_dict.values() if value)

    member_types = list(filter(lambda x: x is not None, [at1, at2, at3, at4]))
    member_types = ['*' if mem_type == '' else mem_type for mem_type in member_types]
    for member in member_types:
        if member == '*':
            continue
        if member not in ref_dict and member not in atom_classes_set:
            raise ForceFieldParseError('AtomType/AtomClass {} not present in AtomTypes reference in the xml'.format(member))
    return member_types


def _parse_default_units(unit_tag):
    if unit_tag is None:
        unit_tag = {}
    units_map = {
        'energy': u.kcal / u.mol,
        'distance': u.nm,
        'mass': u.gram / u.mol,
        'charge': u.coulomb,
        'time': u.ps,
        'temperature': u.K,
        'angle': u.rad
    }
    for attrib, val in unit_tag.items():
        units_map[attrib] = _parse_unit_string(val)
    return units_map


def validate(xml_path_or_etree, schema=None):
    """Validate a given xml file or etree.ElementTree with a reference schema"""
    if schema is None:
        schema_path = os.path.join(os.path.split(os.path.abspath(__file__))[0], 'schema', 'ff-gmso.xsd')
        xml_doc = etree.parse(schema_path)
        xml_schema = etree.XMLSchema(xml_doc)
    else:
        xml_schema = schema

    ff_xml = xml_path_or_etree
    if not isinstance(xml_path_or_etree, etree._ElementTree):
        ff_xml = etree.parse(xml_path_or_etree)

    xml_schema.assertValid(ff_xml)


def _parse_scaling_factors(meta_tag):
    """Parse the scaling factors from the schema"""
    assert meta_tag.tag == 'FFMetaData', 'Can only parse metadata from FFMetaData tag'
    scaling_factors = {'electrostatics14Scale': meta_tag.get('electrostatics14Scale', 1.0),
                       'nonBonded14Scale': meta_tag.get('nonBonded14Scale', 1.0)}
    for key in scaling_factors:
        if type(scaling_factors[key]) != float:
            scaling_factors[key] = float(scaling_factors[key])
    return scaling_factors


def parse_ff_metadata(element):
    """Parse the metadata (units, quantities etc...) from the forcefield XML"""
    metatypes = ['Units']
    parsers = {
        'Units': _parse_default_units,
        'ScalingFactors': _parse_scaling_factors
    }
    ff_meta = {'scaling_factors': parsers['ScalingFactors'](element)}
    for metatype in element:
        if metatype.tag in metatypes:
            ff_meta[metatype.tag] = parsers[metatype.tag](metatype)
    return ff_meta


def parse_ff_atomtypes(atomtypes_el, ff_meta):
    """Given an xml element tree rooted at AtomType, traverse the tree to form a proper topology.core.AtomType"""
    atomtypes_dict = {}
    atomclasses_dict = {}
    ref_dict = ChainMap(atomtypes_dict, atomclasses_dict)
    units_dict = ff_meta['Units']
    atom_types_expression = atomtypes_el.attrib.get('expression', None)
    param_unit_dict = _parse_param_units(atomtypes_el)

    # Parse all the atomTypes and create a new AtomType
    for atom_type in atomtypes_el.getiterator('AtomType'):
        ctor_kwargs = {
            'name': 'AtomType',
            'mass': 0.0 * u.g / u.mol,
            'expression': '4*epsilon*((sigma/r)**12 - (sigma/r)**6)',
            'parameters': None,
            'charge': 0.0 * u.elementary_charge,
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
        if isinstance(ctor_kwargs['charge'], str):
            ctor_kwargs['charge'] = u.unyt_quantity(float(ctor_kwargs['charge']), units_dict['charge'])
        params_dict = _parse_params_values(atom_type, param_unit_dict, 'AtomType')
        if not ctor_kwargs['parameters'] and params_dict:
            ctor_kwargs['parameters'] = params_dict
            valued_param_vars = set(sympify(param) for param in params_dict.keys())
            ctor_kwargs['independent_variables'] = sympify(atom_types_expression).free_symbols - valued_param_vars

        _check_valid_string(ctor_kwargs['name'])
        this_atom_type = AtomType(**ctor_kwargs)
        atomtypes_dict[this_atom_type.name] = this_atom_type
        if this_atom_type.atomclass:
            atomclass_atom_types = atomclasses_dict.get(this_atom_type.atomclass, [])
            atomclass_atom_types.append(this_atom_type)
    return ref_dict


TAG_TO_CLASS_MAP = {
    'BondType': BondType,
    'AngleType': AngleType,
    'DihedralType': DihedralType,
    'ImproperType': ImproperType
}


def parse_ff_connection_types(connectiontypes_el, atomtypes_dict, child_tag='BondType'):
    """Given an XML etree Element rooted at BondTypes, parse the XML to create topology.core.AtomTypes,"""
    connectiontypes_dict = {}
    connectiontype_expression = connectiontypes_el.attrib.get('expression', None)
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

        ctor_kwargs['member_types'] = _check_atomtype_existence(connection_type, atomtypes_dict)
        if not ctor_kwargs['parameters']:
            ctor_kwargs['parameters'] = _parse_params_values(connection_type,
                                                             param_unit_dict,
                                                             child_tag,
                                                             ctor_kwargs['expression'])

        valued_param_vars = set(sympify(param) for param in ctor_kwargs['parameters'].keys())
        ctor_kwargs['independent_variables'] = sympify(connectiontype_expression).free_symbols - valued_param_vars
        this_conn_type_key = DICT_KEY_SEPARATOR.join(ctor_kwargs['member_types'])
        this_conn_type = TAG_TO_CLASS_MAP[child_tag](**ctor_kwargs)
        connectiontypes_dict[this_conn_type_key] = this_conn_type

    return connectiontypes_dict


def _parse_unit_string(string):
    """
    Converts a string with unyt units and physical constants to a taggable unit value
    """
    string = string.replace("deg", "__deg")
    string = string.replace("rad", "__rad")
    expr = sympify(str(string))

    sympy_subs = []
    unyt_subs = []

    for symbol in expr.free_symbols:
        try:
            symbol_unit = _unyt_dictionary[symbol.name.strip('_')]
        except KeyError:
            raise u.exceptions.UnitParseError(
                    "Could not find unit symbol",
                    "'{}' in the provided symbols.".format(symbol.name)
                    )
        if isinstance(symbol_unit, u.Unit):
            sympy_subs.append((symbol.name, symbol_unit.base_value))
            unyt_subs.append((symbol.name, symbol_unit.get_base_equivalent().expr))
        elif isinstance(symbol_unit, u.unyt_quantity):
            sympy_subs.append((symbol.name, float(symbol_unit.in_base().value)))
            unyt_subs.append((symbol.name, symbol_unit.units.get_base_equivalent().expr))

    return u.Unit(float(expr.subs(sympy_subs)) * u.Unit(str(expr.subs(unyt_subs))))
