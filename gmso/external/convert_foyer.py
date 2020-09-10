from warnings import warn

import numpy as np
import unyt as u
import os
import sympy
import pathlib

from lxml import etree
from gmso.utils.io import has_foyer  # Need for foyer
from gmso.exceptions import ForceFieldParseError

if has_foyer:
    import foyer


def from_foyer(foyer_xml, gmso_xml=None):
    """Convert a foyer XML to a gmso XML

    Parameters
    ----------
    foyer_xml : XML file
        An XML file in the foyer format 
   
    Returns
    -------
    gmso_xml : XML file, default=None
        An XML file in the gmso format
    """

    if gmso_xml is None:
        stem_file_name = pathlib.Path(str(foyer_xml)).stem
        suffix = "_gmso.xml"
        gmso_xml = pathlib.Path(stem_file_name + suffix)
    else:
        gmso_xml = pathlib.Path(gmso_xml)

    foyer_xml_tree = etree.parse(foyer_xml)
    f_kwargs = {
        "atom_types": [],
        "non_bonded_forces": [],
        "harmonic_bond_types": [],
        "harmonic_angle_types": [],
        "urey_bradley_angle_types": [],
        "rb_torsion_dihedral_types": [],
        "periodic_torsion_dihedral_types": [],
        "periodic_improper_types": [],
    }

    # Try to load in AtomType section
    # Load in AtomTypes section otherwise
    atom_types_el = foyer_xml_tree.findall("AtomType")
    if len(atom_types_el) == 0:
        atom_types_el = foyer_xml_tree.findall("AtomTypes")
        if len(atom_types_el) == 0:
            raise ForceFieldParseError

    for atom_types in atom_types_el:
        for atom_type in atom_types.getiterator("Type"):
            f_kwargs["atom_types"].append(atom_type)
    nonbonded_force_el = foyer_xml_tree.findall("NonbondedForce")
    for atom_types in nonbonded_force_el:
        f_kwargs["coulomb14scale"] = atom_types.attrib["coulomb14scale"]
        f_kwargs["lj14scale"] = atom_types.attrib["lj14scale"]
        for atom_type in atom_types.getiterator("Atom"):
            f_kwargs["non_bonded_forces"].append(atom_type)

    harmonic_bond_force_el = foyer_xml_tree.findall("HarmonicBondForce")
    for hbf in harmonic_bond_force_el:
        for bond_type in hbf.getiterator("Bond"):
            f_kwargs["harmonic_bond_types"].append(bond_type)

    harmonic_angle_force_el = foyer_xml_tree.findall("HarmonicAngleForce")
    for haf in harmonic_angle_force_el:
        for angle_type in haf.getiterator("Angle"):
            f_kwargs["harmonic_angle_types"].append(angle_type)

    urey_bradley_angle_el = foyer_xml_tree.findall("AmoebaUreyBradleyForce")
    for ubf in urey_bradley_angle_el:
        for angle_type in ubf.getiterator("UreyBradley"):
            f_kwargs["urey_bradley_angle_types"].append(angle_type)

    periodic_torsion_force_el = foyer_xml_tree.findall("PeriodicTorsionForce")
    for ptf in periodic_torsion_force_el:
        for dihedral_type in ptf.getiterator("Proper"):
            f_kwargs["periodic_torsion_dihedral_types"].append(dihedral_type)
        for dihedral_type in ptf.getiterator("Improper"):
            f_kwargs["periodic_improper_types"].append(dihedral_type)

    rb_torsion_force_el = foyer_xml_tree.findall("RBTorsionForce")
    for rbf in rb_torsion_force_el:
        for dihedral_type in rbf.getiterator("Proper"):
            f_kwargs["rb_torsion_dihedral_types"].append(dihedral_type)

    _write_gmso_xml(gmso_xml, **f_kwargs)


def _write_gmso_xml(gmso_xml, **kwargs):
    """Given the set of keyword arguments, write a gmso topology's forcefield xml file"""
    ff_kwargs = {
        "atom_types": kwargs.get("atom_types", []),
        "non_bonded_forces": kwargs.get("non_bonded_forces", []),
        "harmonic_bond_types": kwargs.get("harmonic_bond_types", []),
        "harmonic_angle_types": kwargs.get("harmonic_angle_types", []),
        "urey_bradley_angle_types": kwargs.get("urey_bradley_angle_types", []),
        "periodic_torsion_dihedral_types": kwargs.get(
            "periodic_torsion_dihedral_types", []
        ),
        "periodic_improper_types": kwargs.get("periodic_improper_types", []),
        "rb_torsion_dihedral_types": kwargs.get("rb_torsion_dihedral_types", []),
    }
    forceField = etree.Element("ForceField")
    forceField.attrib["name"] = pathlib.Path(str(gmso_xml)).stem
    forceField.attrib["version"] = "0.0.1"

    ffMeta = _create_subelement(forceField, "FFMetaData")
    units = _create_subelement(
        ffMeta,
        name="Units",
        attrib_dict={
            "energy": "kJ/mol",
            "mass": "amu",
            "charge": "coulomb",
            "distance": "nm",
        },
    )
    # AtomTypes
    nonBondedAtomTypes = _create_subelement(
        forceField,
        "AtomTypes",
        attrib_dict={
            "expression": "ep * ((sigma/r)**12 - (sigma/r)**6) + q / (e0 * r)",
        },
    )

    # NonBondedForces
    nonBondedAtomTypesParamsUnitsDef_ep = _create_subelement(
        nonBondedAtomTypes,
        "ParametersUnitDef",
        attrib_dict={"parameter": "ep", "unit": "kJ/mol",},
    )
    nonBondedAtomTypesParamsUnitsDef_sigma = _create_subelement(
        nonBondedAtomTypes,
        "ParametersUnitDef",
        attrib_dict={"parameter": "sigma", "unit": "nm",},
    )
    nonBondedAtomTypesParamsUnitsDef_e0 = _create_subelement(
        nonBondedAtomTypes,
        "ParametersUnitDef",
        attrib_dict={"parameter": "e0", "unit": "A**2*s**4/(kg*m**3)",},
    )
    nonBondedAtomTypesParamsUnitsDef_q = _create_subelement(
        nonBondedAtomTypes,
        "ParametersUnitDef",
        attrib_dict={"parameter": "q", "unit": "coulomb",},
    )

    for atom_type in ff_kwargs["atom_types"]:
        thisAtomType = _create_subelement(
            nonBondedAtomTypes,
            "AtomType",
            attrib_dict={
                "name": atom_type.get("name", "AtomType"),
                "atomclass": atom_type.get("class", ""),
                "element": atom_type.get("element", ""),
                "charge": atom_type.get("charge", "0.0"),
                "mass": atom_type.get("mass", "0.0"),
                "definition": atom_type.get("def", ""),
                "description": atom_type.get("desc", ""),
            },
        )

    for i, atom_type in enumerate(ff_kwargs["non_bonded_forces"]):
        thisAtomType = nonBondedAtomTypes.find(
            './/AtomType[@name="{}"]'.format(atom_type.get("type"))
        )
        thisAtomType.attrib["name"] = atom_type.get("type", "AtomType")
        parameters = etree.SubElement(thisAtomType, "Parameters")
        parameter_ep = _create_subelement(
            parameters,
            "Parameter",
            attrib_dict={"name": "ep", "value": atom_type.get("epsilon"),},
        )
        parameter_sigma = _create_subelement(
            parameters,
            "Parameter",
            attrib_dict={"name": "sigma", "value": atom_type.get("sigma"),},
        )
        parameter_e0 = _create_subelement(
            parameters, "Parameter", attrib_dict={"name": "e0", "value": "8.8542e-12",}
        )
        parameter_q = _create_subelement(
            parameters,
            "Parameter",
            attrib_dict={"name": "q", "value": atom_type.get("charge"),},
        )

    # HarmonicBondTypes
    harmonicBondTypes = _create_subelement(
        forceField, "BondTypes", attrib_dict={"expression": "k * (r-r_eq)**2",}
    )

    harmonicBondTypesParamsUnitsDef_k = _create_subelement(
        harmonicBondTypes,
        "ParametersUnitDef",
        attrib_dict={"parameter": "k", "unit": "kJ/nm**2",},
    )

    harmonicBondTypesParamsUnitsDef_r_eq = _create_subelement(
        harmonicBondTypes,
        "ParametersUnitDef",
        attrib_dict={"parameter": "r_eq", "unit": "nm"},
    )

    # HarmonicAngleTypes
    harmonicAngleTypes = _create_subelement(
        forceField,
        "AngleTypes",
        attrib_dict={"expression": "k * (theta - theta_eq)**2",},
    )

    harmonicAngleTypesParamsUnitsDef_k = _create_subelement(
        harmonicAngleTypes,
        "ParametersUnitDef",
        attrib_dict={"parameter": "k", "unit": "kJ/radian**2",},
    )

    harmonicBondTypesParamsUnitsDef_theta_eq = _create_subelement(
        harmonicAngleTypes,
        "ParametersUnitDef",
        attrib_dict={"parameter": "theta_eq", "unit": "radian",},
    )

    for i, bond_type in enumerate(ff_kwargs["harmonic_bond_types"]):
        thisBondType = _create_subelement(
            harmonicBondTypes,
            "BondType",
            attrib_dict={
                "name": bond_type.get("name", "BondType-Harmonic-{}".format(i + 1)),
            },
        )
        for j, item in enumerate(bond_type.items()):
            if "type" in item[0]:
                thisBondType.attrib["type{}".format(j + 1)] = bond_type.get(
                    "type{}".format(j + 1), "c{}".format(j + 1)
                )
            elif "class" in item[0]:
                thisBondType.attrib["type{}".format(j + 1)] = bond_type.get(
                    "class{}".format(j + 1), "c{}".format(j + 1)
                )

        parameters = etree.SubElement(thisBondType, "Parameters")
        parameters_k = _create_subelement(
            parameters,
            "Parameter",
            attrib_dict={"name": "k", "value": bond_type.get("k", "1.0"),},
        )

        parameters_k = _create_subelement(
            parameters,
            "Parameter",
            attrib_dict={"name": "r_eq", "value": bond_type.get("length", "1.0"),},
        )

    # UreyBradleyAngleTypes
    ureybradleyAngleTypes = _create_subelement(
        forceField, "AngleTypes", attrib_dict={"expression": "k * (w - w_0) ** 2",}
    )

    ureybradleyAngleTypesParamsUnitsDef_k = _create_subelement(
        ureybradleyAngleTypes,
        "ParametersUnitDef",
        attrib_dict={"parameter": "k", "unit": "kJ/radian**2"},
    )

    ureybradleyAngleTypesParamsUnitsDef_w_0 = _create_subelement(
        ureybradleyAngleTypes,
        "ParametersUnitDef",
        attrib_dict={"parameter": "w_0", "unit": "nm"},
    )

    for i, angle_type in enumerate(ff_kwargs["harmonic_angle_types"]):
        thisAngleType = _create_subelement(
            harmonicAngleTypes,
            "AngleType",
            attrib_dict={
                "name": angle_type.get("name", "AngleType-Harmonic-{}".format(i + 1)),
            },
        )
        for j, item in enumerate(angle_type.items()):
            if "type" in item[0]:
                thisAngleType.attrib["type{}".format(j + 1)] = angle_type.get(
                    "type{}".format(j + 1), "c{}".format(j + 1)
                )
            elif "class" in item[0]:
                thisAngleType.attrib["type{}".format(j + 1)] = angle_type.get(
                    "class{}".format(j + 1), "c{}".format(j + 1)
                )

        parameters = etree.SubElement(thisAngleType, "Parameters")
        parameter_k = _create_subelement(
            parameters,
            "Parameter",
            attrib_dict={"name": "k", "value": angle_type.get("k", "1.0"),},
        )

        parameter_r_eq = _create_subelement(
            parameters,
            "Parameter",
            attrib_dict={"name": "theta_eq", "value": angle_type.get("angle", "1.0"),},
        )

    for i, angle_type in enumerate(ff_kwargs["urey_bradley_angle_types"]):
        thisAngleType = _create_subelement(
            ureybradleyAngleTypes,
            "AngleType",
            attrib_dict={
                "name": angle_type.get(
                    "name", "AngleType-UreyBradley-{}".format(i + 1)
                ),
                "type1": angle_type.get("type1", "t1"),
                "type2": angle_type.get("type2", "t2"),
                "type3": angle_type.get("type3", "t3"),
            },
        )

        parameters = etree.SubElement(thisAngleType, "Parameters")
        parameter_k = _create_subelement(
            parameters,
            "Parameter",
            attrib_dict={"name": "k", "value": angle_type.get("k", "1.0"),},
        )
        parameter_r_eq = _create_subelement(
            parameters,
            "Parameter",
            attrib_dict={"name": "w_0", "value": angle_type.get("d", "1.0"),},
        )

    # PeriodicTorsionDihedralTypes and PeriodicImproperDihedralTypes
    if len(ff_kwargs['periodic_torsion_dihedral_types']) > 0:
        periodicTorsionDihedralTypes = _create_subelement(
            forceField,
            "DihedralTypes",
            attrib_dict={"expression": "k * (1 + cos(n * phi - delta))",},
        )

        max_j = 0
        for i, dihedral_type in enumerate(ff_kwargs["periodic_torsion_dihedral_types"]):
            thisDihedralType = _create_subelement(
                periodicTorsionDihedralTypes,
                "DihedralType",
                attrib_dict={
                    "name": dihedral_type.get(
                        "name", "DihedralType-Periodic-Proper-{}".format(i + 1)
                    ),
                },
            )

            for j, item in enumerate(dihedral_type.items()):
                if "type" in item[0]:
                    thisDihedralType.attrib["type{}".format(j + 1)] = dihedral_type.get(
                        "type{}".format(j + 1), "c{}".format(j + 1)
                    )
                elif "class" in item[0]:
                    thisDihedralType.attrib["type{}".format(j + 1)] = dihedral_type.get(
                        "class{}".format(j + 1), "c{}".format(j + 1)
                    )

            parameters = etree.SubElement(thisDihedralType, "Parameters")

            j = 1
            while dihedral_type.get("k{}".format(j)):
                parameter_k_name = "k{}".format(j)
                parameter_k = _create_subelement(
                    parameters,
                    "Parameter",
                    attrib_dict={
                        "name": parameter_k_name,
                        "value": dihedral_type.get(parameter_k_name),
                    },
                )
                parameter_n = _create_subelement(
                    parameters,
                    "Parameter",
                    attrib_dict={
                        "name": "n{}".format(j),
                        "value": dihedral_type.get("periodicity{}".format(j)),
                    },
                )

                parameter_delta = _create_subelement(
                    parameters,
                    "Parameter",
                    attrib_dict={
                        "name": "delta{}".format(j),
                        "value": dihedral_type.get("phase{}".format(j)),
                    },
                )
                j += 1

            if j > max_j:
                max_j = j
        for k in range(0, max_j):
            periodicTorsionDihedralTypesParamsUnitsDef_k = etree.Element(
                "ParametersUnitDef"
            )
            periodicTorsionDihedralTypesParamsUnitsDef_k.attrib["parameter"] = "k{}".format(
                k
            )
            periodicTorsionDihedralTypesParamsUnitsDef_k.attrib["unit"] = "kJ"
            periodicTorsionDihedralTypes.insert(
                0, periodicTorsionDihedralTypesParamsUnitsDef_k
            )

            periodicTorsionDihedralTypesParamsUnitsDef_n = etree.Element(
                "ParametersUnitDef"
            )
            periodicTorsionDihedralTypesParamsUnitsDef_n.attrib["parameter"] = "n{}".format(
                k
            )
            periodicTorsionDihedralTypesParamsUnitsDef_n.attrib["unit"] = "dimensionless"
            periodicTorsionDihedralTypes.insert(
                0, periodicTorsionDihedralTypesParamsUnitsDef_n
            )

            periodicTorsionDihedralTypesParamsUnitsDef_del = etree.Element(
                "ParametersUnitDef"
            )
            periodicTorsionDihedralTypesParamsUnitsDef_del.attrib[
                "parameter"
            ] = "delta{}".format(k)
            periodicTorsionDihedralTypesParamsUnitsDef_del.attrib["unit"] = "radian"
            periodicTorsionDihedralTypes.insert(
                0, periodicTorsionDihedralTypesParamsUnitsDef_del
            )

    elif len(ff_kwargs['periodic_improper_types']) > 0:
        max_j = 0
        periodicImproperTypes = _create_subelement(
            forceField,
            "DihedralTypes",
            attrib_dict={"expression": "k * (1 + cos(n * phi - delta))",},
        )
        for i, dihedral_type in enumerate(ff_kwargs["periodic_improper_types"]):
            thisDihedralType = _create_subelement(
                periodicImproperTypes,
                "DihedralType",
                attrib_dict={
                    "name": dihedral_type.get(
                        "name", "DihedralType-Periodic-Improper-{}".format(i + 1)
                    ),
                },
            )

            for j, item in enumerate(dihedral_type.items()):
                if "type" in item[0]:
                    thisDihedralType.attrib["type{}".format(j + 1)] = dihedral_type.get(
                        "type{}".format(j + 1), "c{}".format(j + 1)
                    )
                elif "class" in item[0]:
                    thisDihedralType.attrib["type{}".format(j + 1)] = dihedral_type.get(
                        "class{}".format(j + 1), "c{}".format(j + 1)
                    )

            parameters = etree.SubElement(thisDihedralType, "Parameters")

            j = 1
            while dihedral_type.get("k{}".format(j)):
                parameter_k_name = "k{}".format(j)
                parameter_k = _create_subelement(
                    parameters,
                    "Parameter",
                    attrib_dict={
                        "name": parameter_k_name,
                        "value": dihedral_type.get(parameter_k_name),
                    },
                )
                parameter_n = _create_subelement(
                    parameters,
                    "Parameter",
                    attrib_dict={
                        "name": "n{}".format(j),
                        "value": dihedral_type.get("periodicity{}".format(j)),
                    },
                )
                parameter_delta = _create_subelement(
                    parameters,
                    "Parameter",
                    attrib_dict={
                        "name": "delta{}".format(j),
                        "value": dihedral_type.get("phase{}".format(j)),
                    },
                )
                j += 1
            if j > max_j:
                max_j = j

        for k in range(0, max_j):
            periodicImproperTypesParamsUnitsDef_k = etree.Element("ParametersUnitDef")
            periodicImproperTypesParamsUnitsDef_k.attrib["parameter"] = "k{}".format(k)
            periodicImproperTypesParamsUnitsDef_k.attrib["unit"] = "kJ"
            periodicImproperTypes.insert(0, periodicImproperTypesParamsUnitsDef_k)

            periodicImproperTypesParamsUnitsDef_n = etree.Element("ParametersUnitDef")
            periodicImproperTypesParamsUnitsDef_n.attrib["parameter"] = "n{}".format(k)
            periodicImproperTypesParamsUnitsDef_n.attrib["unit"] = "dimensionless"
            periodicImproperTypes.insert(0, periodicImproperTypesParamsUnitsDef_n)

            periodicImproperTypesParamsUnitsDef_del = etree.Element("ParametersUnitDef")
            periodicImproperTypesParamsUnitsDef_del.attrib["parameter"] = "delta{}".format(
                k
            )
            periodicImproperTypesParamsUnitsDef_del.attrib["unit"] = "degree"
            periodicImproperTypes.insert(0, periodicImproperTypesParamsUnitsDef_del)

    # RBTorsionDihedralTypes
    rbTorsionDihedralTypes = _create_subelement(
        forceField,
        "DihedralTypes",
        attrib_dict={
            "expression": "c0 * cos(phi)**0 + c1 * cos(phi)**1 + c2 * cos(phi)**2 + c3 * cos(phi)**3 + c4 * cos(phi)**4 + c5 * cos(phi)**5",
        },
    )

    max_j = 0
    for i, dihedral_type in enumerate(ff_kwargs["rb_torsion_dihedral_types"]):
        thisDihedralType = _create_subelement(
            rbTorsionDihedralTypes,
            "DihedralType",
            attrib_dict={
                "name": dihedral_type.get(
                    "name", "DihedralType-RB-Proper-{}".format(i + 1)
                ),
            },
        )

        for j, item in enumerate(dihedral_type.items()):
            if "type" in item[0]:
                thisDihedralType.attrib["type{}".format(j + 1)] = dihedral_type.get(
                    "type{}".format(j + 1), "c{}".format(j + 1)
                )
            elif "class" in item[0]:
                thisDihedralType.attrib["type{}".format(j + 1)] = dihedral_type.get(
                    "class{}".format(j + 1), "c{}".format(j + 1)
                )

        parameters = etree.SubElement(thisDihedralType, "Parameters")

        j = 0
        while dihedral_type.get("c{}".format(j)):
            parameter_c_name = "c{}".format(j)
            parameter_c = _create_subelement(
                parameters,
                "Parameter",
                attrib_dict={
                    "name": parameter_c_name,
                    "value": dihedral_type.get(parameter_c_name),
                },
            )
            j += 1
        if j > max_j:
            max_j = j
    for k in range(0, max_j):
        rbTorsionDihedralTypesParamsUnitsDef_c = etree.Element("ParametersUnitDef")
        rbTorsionDihedralTypesParamsUnitsDef_c.attrib["parameter"] = "c{}".format(k)
        rbTorsionDihedralTypesParamsUnitsDef_c.attrib["unit"] = "kJ/mol"
        rbTorsionDihedralTypes.insert(0, rbTorsionDihedralTypesParamsUnitsDef_c)

    ff_tree = etree.ElementTree(forceField)
    ff_tree.write(
        str(gmso_xml), pretty_print=True, xml_declaration=True, encoding="utf-8"
    )


def _create_subelement(root_el, name, attrib_dict=None):
    sub_el = etree.SubElement(root_el, name, attrib_dict)
    return sub_el
