import os

import unyt as u
from sympy import sympify
from lxml import etree

from topology.core.atom_type import AtomType
from topology.exceptions import ForceFieldParseError

__all__ = ['validate', 'parse_ff_metadata', 'parse_ff_atomtype', ]


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
        schema_path = os.path.join(os.path.split(os.path.abspath(__file__))[0], 'schema', 'article-schema.xsd')

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


def parse_ff_atomtype(atomtypes_el, meta_map=None):
    """Given an xml element tree rooted at AtomType, traverse the tree to form a proper topology.core.AtomType"""
    # First of all Insert ParamUnits into a default dictionary if exists
    atomtypes_dict = {}
    param_unit_dict = {}
    units_dict = meta_map['Units']
    atom_types_expression = atomtypes_el.attrib['expression']
    for param_unit in atomtypes_el.getiterator('ParametersUnitDef'):
        param_unit_dict[param_unit.attrib['parameter']] = param_unit.attrib['unit']
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

    for atom_type in atomtypes_el.getiterator('AtomType'):
        for kwarg in ctor_kwargs.keys():
            ctor_kwargs[kwarg] = atom_type.attrib.get(kwarg, ctor_kwargs[kwarg])
        if isinstance(ctor_kwargs['mass'], str):
            ctor_kwargs['mass'] = u.unyt_quantity(float(ctor_kwargs['mass']), units_dict['mass'])
        if isinstance(ctor_kwargs['overrides'], str):
            ctor_kwargs['overrides'] = set(ctor_kwargs['overrides'].split(','))

        # Tag of type Parameters can exist atmost once
        params_dict = {}
        for param in atom_type.find('Parameters').getiterator('Parameter'):
            if param.attrib['name'] not in param_unit_dict:
                raise ForceFieldParseError('Parameters {} with Unknown units found'.format(param.attrib['name']))
            param_name = param.attrib['name']
            param_unit = param_unit_dict[param_name]
            param_value = u.unyt_quantity(float(param.attrib['value']), param_unit)
            params_dict[param_name] = param_value

        if not ctor_kwargs['parameters']:
            ctor_kwargs['parameters'] = params_dict
        ctor_kwargs['independent_variables'] = sympify(atom_types_expression).free_symbols - set(params_dict.keys())
        this_atom_type = AtomType(**ctor_kwargs)
        atomtypes_dict[this_atom_type] = this_atom_type
        print(atomtypes_dict)
    return atomtypes_dict

