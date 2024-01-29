"""Utility module for operating on gmso.ForceFields."""

import os
import re

import numpy as np
import unyt as u
from lxml import etree
from sympy import sympify

from gmso.core.angle_type import AngleType
from gmso.core.atom_type import AtomType
from gmso.core.bond_type import BondType
from gmso.core.dihedral_type import DihedralType
from gmso.core.improper_type import ImproperType
from gmso.core.pairpotential_type import PairPotentialType
from gmso.exceptions import (
    ForceFieldError,
    ForceFieldParseError,
    MissingAtomTypesError,
    MixedClassAndTypesError,
)
from gmso.utils._constants import FF_TOKENS_SEPARATOR

__all__ = [
    "validate",
    "parse_ff_metadata",
    "parse_ff_atomtypes",
    "parse_ff_connection_types",
    "parse_ff_pairpotential_types",
]


# Create a dictionary of units
_unyt_dictionary = {}
for name, item in vars(u).items():
    if isinstance(item, u.Unit) or isinstance(item, u.unyt_quantity):
        _unyt_dictionary.update({name: item})


def _check_valid_string(type_str):
    """Verify that the strings in the forcefield are valid."""
    if FF_TOKENS_SEPARATOR in type_str:
        raise ForceFieldError(
            "Please do not use {} in type string".format(FF_TOKENS_SEPARATOR)
        )


def _parse_param_units(parent_tag):
    """Return the units associated with parameters."""
    param_unit_dict = {}
    params_iter = parent_tag.getiterator("ParametersUnitDef")
    for param_unit in params_iter:
        param_unit_dict[param_unit.attrib["parameter"]] = _parse_unit_string(
            param_unit.attrib["unit"]
        )
    return param_unit_dict


def _parse_params_values(parent_tag, units_dict, child_tag, expression=None):
    """Return the values of the parameters."""
    # Tag of type Parameters can exist atmost once
    params_dict = {}
    if parent_tag.find("Parameters") is None:
        return params_dict
    for param in parent_tag.find("Parameters").getiterator("Parameter"):
        if param.attrib["name"] not in units_dict:
            raise ForceFieldParseError(
                "Parameters {} with Unknown units found".format(
                    param.attrib["name"]
                )
            )
        param_name = param.attrib["name"]
        param_unit = units_dict[param_name]
        if param.attrib.get("value"):
            param_value = u.unyt_quantity(
                float(param.attrib["value"]), param_unit
            )
        else:
            children = param.getchildren()
            if len(children) == 0:
                raise ForceFieldParseError(
                    f"Neither a single value nor a sequence of values "
                    f"is specified for parameter {param_name}, please specify "
                    f"either a single value as an attribute value or a sequence "
                    f"of values."
                )
            value_array = np.array(
                [value.text for value in children], dtype=float
            )
            param_value = u.unyt_array(value_array, param_unit)

        params_dict[param_name] = param_value
    param_ref_dict = units_dict
    if child_tag == "DihedralType":
        if not expression:
            raise ForceFieldError(
                "Cannot consolidate parameters without an expression"
            )
        _consolidate_params(params_dict, expression)
        param_ref_dict = _consolidate_params(
            units_dict, expression, update_orig=False
        )

    for param in param_ref_dict:
        if param not in params_dict:
            raise ForceFieldParseError(
                "Parameter {} is in units but cannot be found in parameters list".format(
                    param
                )
            )
    return params_dict


def _consolidate_params(params_dict, expression, update_orig=True):
    """Consolidate parameters for a specific expression."""
    to_del = []
    new_dict = {}
    match_string = "|".join(
        str(symbol) for symbol in sympify(expression).free_symbols
    )
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


def _get_member_types(tag):
    """Return the types of the members, handle wildcards."""
    at1 = tag.attrib.get("type1")
    at2 = tag.attrib.get("type2")
    at3 = tag.attrib.get("type3")
    at4 = tag.attrib.get("type4")

    member_types = filter(lambda x: x is not None, [at1, at2, at3, at4])
    member_types = [
        "*" if mem_type == "" else mem_type for mem_type in member_types
    ]
    return member_types or None


def _get_member_classes(tag):
    """Return the classes of the members, handle wildcards."""
    at1 = tag.attrib.get("class1")
    at2 = tag.attrib.get("class2")
    at3 = tag.attrib.get("class3")
    at4 = tag.attrib.get("class4")

    member_classes = filter(lambda x: x is not None, [at1, at2, at3, at4])
    member_classes = [
        "*" if mem_type == "" else mem_type for mem_type in member_classes
    ]

    return member_classes or None


def _parse_default_units(unit_tag):
    """Check if units are different from the default set below."""
    if unit_tag is None:
        unit_tag = {}
    units_map = {
        "energy": u.kcal / u.mol,
        "distance": u.nm,
        "mass": u.gram / u.mol,
        "charge": u.coulomb,
        "time": u.ps,
        "temperature": u.K,
        "angle": u.rad,
    }
    for attrib, val in unit_tag.items():
        units_map[attrib] = _parse_unit_string(val)
    return units_map


def validate(gmso_xml_or_etree, strict=True, greedy=True):
    """Validate the gmso XML file or etree.ElementTree.

    This function validates the given gmso XML file or etree.ElementTree object
    against the gmso XML schema and optionally provide additional validation in
    strict mode.

    Parameters
    ----------
    gmso_xml_or_etree: str, pathlib.Path, etree.ElementTree
        The XML file to perform validation for
    strict: bool, default=True
        If true, perform a strict validation which includes:
            1. Check if all the atom_types/classes for individual entries
             in the `BondTypes`, `AngleTypes`, `DihedralTypes` section are
             found in the entries of the `AtomTypes` section
    greedy: bool, default=False
        If true, report all the mismatches that were found in the strict
        validation else fail in the first mismatch

    Notes
    -----
    `greedy` is only used when `strict` is True.

    See Also
    --------
    _validate_schema
        Validates a xml file or etree.ElementTree with a reference schema

    Raises
    ------
    MissingAtomTypesError
        If `strict` is True and all the atom_types/classes for individual entries
        in the `BondTypes`, `AngleTypes`, `DihedralTypes` section are found in
        the entries of the `AtomTypes` section
    """
    ff_etree = _validate_schema(xml_path_or_etree=gmso_xml_or_etree)
    _assert_membertype_class_exclusivity(ff_etree)
    if strict:
        missing = _find_missing_atom_types_or_classes(ff_etree, greedy=greedy)
        if missing:
            raise MissingAtomTypesError(
                f"Atom types/classes {missing} are missing in the AtomTypes "
                f"section but present in the BondTypes/AngleTypes/DihedralTypes "
                f"section of the ForceField XML file. If this behavior is intended, "
                f"please disable this check by setting strict=False."
            )


def _assert_membertype_class_exclusivity(root):
    """Check if there's a criss-cross between type and class in Bonded Potentials."""
    for idx, potential_tag in enumerate(
        ["BondType", "AngleType", "DihedralType", "ImproperType"], start=2
    ):
        potential_iter = root.iterfind(f".//{potential_tag}")
        for potential in potential_iter:
            if potential_tag == "ImproperType":
                iter_end_idx = idx
            else:
                iter_end_idx = idx + 1
            types_and_classes = (
                (
                    potential.attrib.get(f"type{j}"),
                    potential.attrib.get(f"class{j}"),
                )
                for j in range(1, iter_end_idx)
            )
            error_msg = (
                f"{potential_tag} {potential.attrib['name']} has a mix "
                f"of atom type and atom classes "
                f"which is not allowed, please use uniform attributes.\n"
                f"{etree.tostring(potential, encoding='utf-8', pretty_print=True).decode()}"
            )

            types, classes = zip(*types_and_classes)
            if any(types) and any(classes):
                raise MixedClassAndTypesError(error_msg)


def _find_missing_atom_types_or_classes(ff_etree, greedy=False):
    """Iterate through the forcefield tree and find any missing atomtypes or classes."""
    atom_types_iter = ff_etree.iterfind(".//AtomType")
    atom_types_and_classes = set()
    for atom_type in atom_types_iter:
        if atom_type.attrib.get("name"):
            atom_types_and_classes.add(atom_type.attrib["name"])
        if atom_type.attrib.get("atomclass"):
            atom_types_and_classes.add(atom_type.attrib["atomclass"])

    remaining_potentials = [
        ff_etree.iterfind(".//BondType"),
        ff_etree.iterfind(".//BondType"),
        ff_etree.iterfind(".//AngleType"),
        ff_etree.iterfind(".//DihedralType"),
    ]
    member_types_or_classes = set()

    # ToDo: This should be made a wildcard, stored globally
    if "*" not in atom_types_and_classes:
        atom_types_and_classes.add("*")

    for potentials_type in remaining_potentials:
        for potential_type in potentials_type:
            types_or_classes = _get_member_types(potential_type)
            if not types_or_classes:
                types_or_classes = _get_member_classes(potential_type)
            for type_or_class in types_or_classes:
                member_types_or_classes.add(type_or_class)

    missing = []
    for type_or_class in member_types_or_classes:
        if type_or_class not in atom_types_and_classes:
            missing.append(type_or_class)
            if missing and greedy:
                break

    return missing


def _validate_schema(xml_path_or_etree, schema=None):
    """Validate a given xml file or etree.ElementTree with a reference schema."""
    if schema is None:
        schema_path = os.path.join(
            os.path.split(os.path.abspath(__file__))[0], "schema", "ff-gmso.xsd"
        )
        xml_doc = etree.parse(schema_path)
        xml_schema = etree.XMLSchema(xml_doc)
    else:
        xml_schema = schema

    ff_xml = xml_path_or_etree
    if not isinstance(xml_path_or_etree, etree._ElementTree):
        ff_xml = etree.parse(xml_path_or_etree)
    try:
        xml_schema.assertValid(ff_xml)
    except etree.DocumentInvalid as ex:
        message = ex.error_log.last_error.message
        line = ex.error_log.last_error.line
        # rewrite error message for constraint violation
        if ex.error_log.last_error.type_name == "SCHEMAV_CVC_IDC":
            raise ForceFieldParseError(message, line)
    return ff_xml


def _parse_scaling_factors(meta_tag):
    """Parse the scaling factors from the schema."""
    assert (
        meta_tag.tag == "FFMetaData"
    ), "Can only parse metadata from FFMetaData tag"
    scaling_factors = {
        "electrostatics14Scale": meta_tag.get("electrostatics14Scale", 1.0),
        "nonBonded14Scale": meta_tag.get("nonBonded14Scale", 1.0),
    }
    for key in scaling_factors:
        if type(scaling_factors[key]) != float:
            scaling_factors[key] = float(scaling_factors[key])
    return scaling_factors


def parse_ff_metadata(element):
    """Parse the metadata (units, quantities etc...) from the forcefield XML."""
    metatypes = ["Units"]
    parsers = {
        "Units": _parse_default_units,
        "ScalingFactors": _parse_scaling_factors,
    }
    ff_meta = {
        "scaling_factors": parsers["ScalingFactors"](element),
        "combining_rule": element.get("combiningRule", "geometric"),
    }
    for metatype in element:
        if metatype.tag in metatypes:
            ff_meta[metatype.tag] = parsers[metatype.tag](metatype)
    return ff_meta


def parse_ff_atomtypes(atomtypes_el, ff_meta):
    """Given an xml element tree rooted at AtomType, traverse the tree to form a proper topology.core.AtomType."""
    atomtypes_dict = {}
    units_dict = ff_meta["Units"]
    atom_types_expression = atomtypes_el.attrib.get("expression", None)
    param_unit_dict = _parse_param_units(atomtypes_el)

    # Parse all the atomTypes and create a new AtomType
    for atom_type in atomtypes_el.getiterator("AtomType"):
        ctor_kwargs = {
            "name": "AtomType",
            "mass": 0.0 * u.g / u.mol,
            "expression": "4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            "parameters": None,
            "charge": 0.0 * u.elementary_charge,
            "independent_variables": None,
            "atomclass": "",
            "doi": "",
            "overrides": set(),
            "definition": "",
            "description": "",
            "element": "",
        }

        if atom_types_expression:
            ctor_kwargs["expression"] = atom_types_expression

        for kwarg in ctor_kwargs.keys():
            ctor_kwargs[kwarg] = atom_type.attrib.get(kwarg, ctor_kwargs[kwarg])

        tags = {"tags": {"element": ctor_kwargs.pop("element", "")}}

        ctor_kwargs.update(tags)

        if isinstance(ctor_kwargs["mass"], str):
            ctor_kwargs["mass"] = u.unyt_quantity(
                float(ctor_kwargs["mass"]), units_dict["mass"]
            )
        if isinstance(ctor_kwargs["overrides"], str):
            ctor_kwargs["overrides"] = set(
                override.strip()
                for override in ctor_kwargs["overrides"].split(",")
            )
        if isinstance(ctor_kwargs["charge"], str):
            ctor_kwargs["charge"] = u.unyt_quantity(
                float(ctor_kwargs["charge"]), units_dict["charge"]
            )
        params_dict = _parse_params_values(
            atom_type, param_unit_dict, "AtomType"
        )
        if not ctor_kwargs["parameters"] and params_dict:
            ctor_kwargs["parameters"] = params_dict
            valued_param_vars = set(
                sympify(param) for param in params_dict.keys()
            )
            ctor_kwargs["independent_variables"] = (
                sympify(atom_types_expression).free_symbols - valued_param_vars
            )

        _check_valid_string(ctor_kwargs["name"])
        this_atom_type = AtomType(**ctor_kwargs)
        atomtypes_dict[this_atom_type.name] = this_atom_type

    return atomtypes_dict


TAG_TO_CLASS_MAP = {
    "BondType": BondType,
    "AngleType": AngleType,
    "DihedralType": DihedralType,
    "ImproperType": ImproperType,
    "PairPotentialType": PairPotentialType,
}


def parse_ff_connection_types(connectiontypes_el, child_tag="BondType"):
    """Parse an XML etree Element rooted at BondTypes to create topology.core.AtomTypes."""
    connectiontypes_dict = {}
    connectiontype_expression = connectiontypes_el.attrib.get(
        "expression", None
    )
    param_unit_dict = _parse_param_units(connectiontypes_el)

    # Parse all the bondTypes and create a new BondType
    for connection_type in connectiontypes_el.getiterator(child_tag):
        ctor_kwargs = {
            "name": child_tag,
            "expression": "0.5 * k * (r-r_eq)**2",
            "parameters": None,
            "independent_variables": None,
            "member_types": None,
            "member_classes": None,
        }
        if connectiontype_expression:
            ctor_kwargs["expression"] = connectiontype_expression

        for kwarg in ctor_kwargs.keys():
            ctor_kwargs[kwarg] = connection_type.attrib.get(
                kwarg, ctor_kwargs[kwarg]
            )

        ctor_kwargs["member_types"] = _get_member_types(connection_type)
        if not ctor_kwargs["member_types"]:
            ctor_kwargs["member_classes"] = _get_member_classes(connection_type)

        if not ctor_kwargs["parameters"]:
            ctor_kwargs["parameters"] = _parse_params_values(
                connection_type,
                param_unit_dict,
                child_tag,
                ctor_kwargs["expression"],
            )

        valued_param_vars = set(
            sympify(param) for param in ctor_kwargs["parameters"].keys()
        )
        ctor_kwargs["independent_variables"] = (
            sympify(connectiontype_expression).free_symbols - valued_param_vars
        )

        this_conn_type_key = FF_TOKENS_SEPARATOR.join(
            ctor_kwargs.get("member_types") or ctor_kwargs.get("member_classes")
        )

        this_conn_type = TAG_TO_CLASS_MAP[child_tag](**ctor_kwargs)
        connectiontypes_dict[this_conn_type_key] = this_conn_type

    return connectiontypes_dict


def parse_ff_pairpotential_types(pairpotentialtypes_el):
    """Given an XML etree Element rooted at PairPotentialTypes, parse the XML to create topology.core.PairPotentialTypes."""
    pairpotentialtypes_dict = {}
    pairpotentialtype_expression = pairpotentialtypes_el.attrib.get(
        "expression", None
    )
    param_unit_dict = _parse_param_units(pairpotentialtypes_el)

    # Parse all the pairpotentialTypes and create a new PairPotentialType
    for pairpotential_type in pairpotentialtypes_el.getiterator(
        "PairPotentialType"
    ):
        ctor_kwargs = {
            "name": PairPotentialType,
            "expression": "4*epsilon*((sigma/r)**12 - (sigma/r)**6)",
            "parameters": None,
            "independent_variables": None,
            "member_types": None,
        }
        if pairpotentialtype_expression:
            ctor_kwargs["expression"] = pairpotentialtype_expression

        for kwarg in ctor_kwargs.keys():
            ctor_kwargs[kwarg] = pairpotential_type.attrib.get(
                kwarg, ctor_kwargs[kwarg]
            )

        ctor_kwargs["member_types"] = _get_member_types(pairpotential_type)
        if not ctor_kwargs["parameters"]:
            ctor_kwargs["parameters"] = _parse_params_values(
                pairpotential_type,
                param_unit_dict,
                "PairPotentialType",
                ctor_kwargs["expression"],
            )

        valued_param_vars = set(
            sympify(param) for param in ctor_kwargs["parameters"].keys()
        )
        ctor_kwargs["independent_variables"] = (
            sympify(pairpotentialtype_expression).free_symbols
            - valued_param_vars
        )
        this_pp_type_key = FF_TOKENS_SEPARATOR.join(ctor_kwargs["member_types"])
        this_pp_type = TAG_TO_CLASS_MAP["PairPotentialType"](**ctor_kwargs)
        pairpotentialtypes_dict[this_pp_type_key] = this_pp_type

    return pairpotentialtypes_dict


def _parse_unit_string(string):
    """Convert a string with unyt units and physical constants to a taggable unit value."""
    string = string.replace("deg", "__deg")
    string = string.replace("rad", "__rad")
    expr = sympify(str(string))

    sympy_subs = []
    unyt_subs = []

    for symbol in expr.free_symbols:
        try:
            symbol_unit = _unyt_dictionary[symbol.name.strip("_")]
        except KeyError:
            raise u.exceptions.UnitParseError(
                "Could not find unit symbol",
                "'{}' in the provided symbols.".format(symbol.name),
            )
        if isinstance(symbol_unit, u.Unit):
            sympy_subs.append((symbol.name, symbol_unit.base_value))
            unyt_subs.append(
                (symbol.name, symbol_unit.get_base_equivalent().expr)
            )
        elif isinstance(symbol_unit, u.unyt_quantity):
            sympy_subs.append((symbol.name, float(symbol_unit.in_base().value)))
            unyt_subs.append(
                (symbol.name, symbol_unit.units.get_base_equivalent().expr)
            )

    return u.Unit(
        float(expr.subs(sympy_subs)) * u.Unit(str(expr.subs(unyt_subs)))
    )
